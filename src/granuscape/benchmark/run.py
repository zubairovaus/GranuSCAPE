"""Repeated complete processing with explicit first-run and failure accounting."""

from datetime import datetime, timezone
import json
from pathlib import Path
import random
import time
import zipfile

import numpy as np
import pandas as pd

from granuscape import __version__
from granuscape.io.czi import sha256
from granuscape.pipeline.run import run_batch, code_hash
from granuscape.qc.settings import resolve_options
from granuscape.segmentation.backends import MODELS
from .environment import inventory

# These are disjoint top-level costs; named analysis/segmentation substeps are nested.
PRIMARY = [
    "input_hash_seconds",
    "import_seconds",
    "inference_seconds",
    "analysis_seconds",
    "hdf5_export_seconds",
    "qc_export_seconds",
    "csv_export_seconds",
    "work_cleanup_seconds",
    "other_seconds",
    "elapsed_seconds",
]


def summarize(images, output, expected_repeats):
    measured = images.loc[images.phase.eq("measured")].copy()
    complete = measured.loc[measured.status.eq("complete_requires_review")]
    counts = complete.groupby("image_id").repeat.nunique()
    all_counts = measured.groupby("image_id").size()
    accepted = counts.loc[
        (counts == expected_repeats)
        & all_counts.reindex(counts.index).eq(expected_repeats)
    ].index
    primary = complete.loc[complete.image_id.isin(accepted)]
    metrics = [
        c
        for c in images.columns
        if c.endswith("_seconds") and c != "batch_wall_seconds"
    ]
    per_image = (
        primary.groupby("image_id")[metrics].median().reset_index()
        if len(primary)
        else pd.DataFrame(columns=["image_id", *metrics])
    )
    per_image.to_csv(output / "benchmark_per_image.csv", index=False)
    rows = []
    for metric in metrics:
        values = pd.to_numeric(per_image[metric], errors="coerce").dropna()
        if values.empty:
            continue
        rows.append(
            dict(
                metric=metric,
                scope="top_level" if metric in PRIMARY else "nested_component",
                images=len(values),
                raw_successful_timings=int(primary[metric].notna().sum()),
                median_seconds=float(values.median()),
                q25_seconds=float(values.quantile(0.25)),
                q75_seconds=float(values.quantile(0.75)),
                min_seconds=float(values.min()),
                max_seconds=float(values.max()),
            )
        )
    pd.DataFrame(
        rows,
        columns=[
            "metric",
            "scope",
            "images",
            "raw_successful_timings",
            "median_seconds",
            "q25_seconds",
            "q75_seconds",
            "min_seconds",
            "max_seconds",
        ],
    ).to_csv(output / "benchmark_summary.csv", index=False)
    report = {
        "expected_repeats": expected_repeats,
        "distinct_images_attempted": int(measured.image_id.nunique()),
        "images_with_all_repeats_complete": len(accepted),
        "failed_warmup_runs": int(
            (
                (images.phase == "warmup")
                & (images.status != "complete_requires_review")
            ).sum()
        ),
        "failed_measured_runs": int(
            (measured.status != "complete_requires_review").sum()
        ),
        "primary_summary_eligible": len(accepted) == measured.image_id.nunique()
        and len(accepted) > 0,
        "aggregation": "Median across repeated measurements of each image, then distribution across images. Repeats are not independent images.",
        "nested_timers": "Do not add nested segmentation/analysis timers to their parent totals. Missing stages remain missing, not zero.",
    }
    (output / "benchmark_status.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )


def benchmark(config, output=None, repeats=3, warmup=1, seed=17):
    if type(repeats) is not int or repeats < 1 or type(warmup) is not int or warmup < 0:
        raise ValueError("repeats must be >=1 and warmup must be >=0")
    files = sorted(
        p
        for p in Path(config["input"]).rglob("*")
        if p.is_file() and p.suffix.lower() == ".czi"
    )
    if not files:
        raise ValueError("No CZI files selected for benchmark")
    if len({p.stem.casefold() for p in files}) != len(files):
        raise ValueError("Duplicate CZI stems")
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    output = (
        Path(output).resolve()
        if output
        else Path(config["output"]).resolve().parent / ("benchmark_" + stamp)
    )
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Benchmark output must be new/empty")
    output.mkdir(parents=True, exist_ok=True)
    config = dict(config, qc=resolve_options(config.get("qc")))
    environment = inventory(config, output)
    (output / "environment.json").write_text(
        json.dumps(environment, indent=2), encoding="utf-8"
    )
    inputs = [
        dict(
            image_id=p.stem,
            source=str(p.resolve()),
            bytes=p.stat().st_size,
            sha256=sha256(p),
        )
        for p in files
    ]
    pd.DataFrame(inputs).to_csv(output / "benchmark_inputs.csv", index=False)
    protocol = {
        "package_version": __version__,
        "source_code_sha256": code_hash(),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "repeats": repeats,
        "warmup_runs": warmup,
        "warmup_input": files[0].stem if warmup else None,
        "seed": seed,
        "measured_image_count": len(files),
        "configuration": config,
        "models_sha256": {p.name: sha256(p) for p in MODELS.glob("*.ilp")},
        "timing_scope": "Each complete run_batch call; excludes initial parent-package import, hardware inventory and initial input-manifest hashes. Includes CZI processing, fresh backend subprocess startup, inference, measurements, HDF5, PNG/CSV export and configured cleanup.",
        "warmup_scope": "Complete first-image runs, separately retained; conditions caches but does not keep models resident. Every image starts new backend subprocesses.",
        "file_safety": "Original CZI files are read-only. Each repetition writes fresh HDF5 and output files; existing results are not reused or overwritten.",
    }
    (output / "benchmark_protocol.json").write_text(
        json.dumps(protocol, indent=2), encoding="utf-8"
    )
    print(
        f"Benchmark: {len(files)} images × {repeats} repeats; {warmup} separate first-image warm-up run(s). Output: {output}",
        flush=True,
    )
    rng = random.Random(seed)
    image_frames = []
    batches = []
    fatal = False
    plan = [("warmup", i + 1, [files[0]]) for i in range(warmup)]
    for i in range(repeats):
        ordered = files.copy()
        rng.shuffle(ordered)
        plan.append(("measured", i + 1, ordered))
    for phase, repeat, selected in plan:
        folder = output / f"{phase}_{repeat:02d}"
        run_config = dict(config, output=str(folder))
        start = time.perf_counter()
        error = ""
        try:
            result = run_batch(run_config, selected_files=selected)
        except Exception as exc:
            result = 1
            error = str(exc)
        wall = time.perf_counter() - start
        batch = dict(
            phase=phase,
            repeat=repeat,
            wall_seconds=wall,
            images_selected=len(selected),
            exit_code=result,
            error=error,
        )
        report = folder / "run_summary.csv"
        if report.is_file():
            frame = pd.read_csv(report)
            for key, value in [("phase", phase), ("repeat", repeat)]:
                frame[key] = value
            attempted = set(frame.image_id.astype(str))
            missing = [
                dict(
                    image_id=p.stem,
                    status="failed",
                    error="No per-image report produced",
                    phase=phase,
                    repeat=repeat,
                )
                for p in selected
                if p.stem not in attempted
            ]
            if missing:
                frame = pd.concat([frame, pd.DataFrame(missing)], ignore_index=True)
        else:
            frame = pd.DataFrame(
                [
                    dict(
                        image_id=p.stem,
                        status="failed",
                        error=error or "Missing run report",
                        phase=phase,
                        repeat=repeat,
                    )
                    for p in selected
                ]
            )
        for h5_path in sorted((folder / "hdf5").glob("*.h5")):
            import h5py

            with h5py.File(h5_path, "r") as h5:
                provenance = json.loads(h5["provenance/json"].asstr()[()])
            backend_path = output / "backend_versions.json"
            if not backend_path.exists():
                backend_path.write_text(
                    json.dumps(provenance.get("backends", {}), indent=2),
                    encoding="utf-8",
                )
            break
        image_frames.append(frame)
        batch["images_complete"] = int(
            (frame.status == "complete_requires_review").sum()
        )
        batch["throughput_images_per_hour"] = (
            batch["images_complete"] / wall * 3600 if result == 0 else np.nan
        )
        batches.append(batch)
        images = pd.concat(image_frames, ignore_index=True)
        images.to_csv(output / "benchmark_images.csv", index=False)
        pd.DataFrame(batches).to_csv(output / "benchmark_batches.csv", index=False)
        summarize(images, output, repeats)
        if result:
            fatal = True
            if phase == "warmup":
                print(
                    "Warm-up failed; benchmark stopped. Resolve the recorded error before a publication benchmark.",
                    flush=True,
                )
                break
    if not (output / "benchmark_status.json").exists():
        (output / "benchmark_status.json").write_text(
            json.dumps(
                {
                    "primary_summary_eligible": False,
                    "reason": "No measured run completed",
                },
                indent=2,
            )
        )
    selected_outputs = [
        "environment.json",
        "backend_versions.json",
        "benchmark_inputs.csv",
        "benchmark_protocol.json",
        "benchmark_images.csv",
        "benchmark_batches.csv",
        "benchmark_per_image.csv",
        "benchmark_summary.csv",
        "benchmark_status.json",
    ]
    with zipfile.ZipFile(
        output / "BENCHMARK_REPORT.zip", "w", zipfile.ZIP_DEFLATED
    ) as archive:
        for name in selected_outputs:
            path = output / name
            if path.is_file():
                archive.write(path, name)
    print(f'Benchmark report: {output / "BENCHMARK_REPORT.zip"}', flush=True)
    return int(fatal)
