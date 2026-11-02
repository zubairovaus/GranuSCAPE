"""Benchmark accounting and timing must preserve numerical results."""

import json
from pathlib import Path
import pandas as pd
import pytest
from granuscape.benchmark.timing import Timer
from granuscape.benchmark.run import summarize, benchmark


def test_elapsed_time_survives_stage_failure(monkeypatch):
    import granuscape.benchmark.timing as timing

    clock = iter([10.0, 12.5])
    monkeypatch.setattr(timing, "perf_counter", lambda: next(clock))
    timer = Timer()
    with pytest.raises(RuntimeError):
        with timer.span("stage"):
            raise RuntimeError("failed stage")
    assert timer.seconds["stage"] == 2.5


def test_repeats_are_aggregated_within_images_and_failures_remain_visible(tmp_path):
    rows = []
    for image, times, statuses in [
        ("A", [1.0, 3.0], ["complete_requires_review"] * 2),
        ("B", [8.0, 12.0], ["complete_requires_review"] * 2),
        ("C", [20.0, 40.0], ["complete_requires_review", "failed"]),
    ]:
        for repeat, (elapsed, status) in enumerate(zip(times, statuses), 1):
            rows.append(
                dict(
                    image_id=image,
                    repeat=repeat,
                    phase="measured",
                    status=status,
                    elapsed_seconds=elapsed,
                )
            )
    summarize(pd.DataFrame(rows), tmp_path, 2)
    summary = pd.read_csv(tmp_path / "benchmark_summary.csv")
    assert summary.iloc[0].median_seconds == 6.0  # Median of image medians 2 and 10.
    assert summary.iloc[0].images == 2
    report = json.loads((tmp_path / "benchmark_status.json").read_text())
    assert report["failed_measured_runs"] == 1
    assert not report["primary_summary_eligible"]


def test_repetitions_use_fresh_outputs_and_original_inputs(tmp_path, monkeypatch):
    import granuscape.benchmark.run as bench

    data = tmp_path / "data"
    data.mkdir()
    for name in ("A", "B"):
        (data / (name + ".czi")).write_bytes(name.encode())
    calls = []

    def run(config, selected_files=None):
        out = Path(config["output"])
        assert not out.exists()
        out.mkdir()
        calls.append((out.name, [p.name for p in selected_files]))
        pd.DataFrame(
            [
                dict(
                    image_id=p.stem,
                    status="complete_requires_review",
                    elapsed_seconds=1.0,
                    import_seconds=0.1,
                    inference_seconds=0.5,
                    analysis_seconds=0.2,
                    qc_export_seconds=0.1,
                )
                for p in selected_files
            ]
        ).to_csv(out / "run_summary.csv", index=False)
        return 0

    monkeypatch.setattr(bench, "run_batch", run)
    monkeypatch.setattr(
        bench, "inventory", lambda config, output: {"test_fixture": True}
    )
    config = dict(input=str(data), output=str(tmp_path / "unused"))
    assert benchmark(config, tmp_path / "benchmark", repeats=2, warmup=1) == 0
    assert len(calls) == 3 and len({name for name, _ in calls}) == 3
    assert len(calls[0][1]) == 1
    assert all(set(files) == {"A.czi", "B.czi"} for _, files in calls[1:])
    assert (data / "A.czi").read_bytes() == b"A"
    assert (tmp_path / "benchmark/BENCHMARK_REPORT.zip").is_file()
    status = json.loads((tmp_path / "benchmark/benchmark_status.json").read_text())
    assert status["primary_summary_eligible"]
