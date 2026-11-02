import importlib
import pkgutil
import h5py
import numpy as np
import pandas as pd
from granuscape.cells.refinement import classify_cells, measure_cells, refine_cells
from granuscape.granules.objects import object_table
from granuscape.io.czi import parse_metadata
from granuscape.io.hdf5 import write_table, read_table
from granuscape.measurements.profiles import cell_metrics
from granuscape.layers.interface import validate_open_curve
from granuscape.pipeline.analysis import analyze_arrays


def fixture():
    y, x = np.indices((256, 256))
    soma = np.zeros((256, 256), np.uint32)
    nuclei = np.zeros_like(soma)
    for sid, cy in enumerate([35, 80, 125, 170, 215], 1):
        soma[(x - 116) ** 2 + (y - cy) ** 2 < 12**2] = sid
        nuclei[(x - 116) ** 2 + (y - cy) ** 2 < 5**2] = sid
    nid = 6
    for cy in range(10, 250, 12):
        for cx in range(145, 250, 12):
            nuclei[(x - cx) ** 2 + (y - cy) ** 2 < 3**2] = nid
            nid += 1
    af = np.full(soma.shape, 20, np.uint16)
    af[soma > 0] = 100
    af[nuclei > 0] = 10
    raw = dict(
        af488=af, dapi=(nuclei > 0).astype(np.uint16) * 200, rhoda=np.full_like(af, 20)
    )
    granules = np.zeros_like(soma)
    granules[30:33, 118:121] = 1
    masks = dict(soma=soma, nuclei=nuclei, g3bp=granules, eif3=granules.copy())
    return raw, masks


def test_numeric_modules_import():
    for name in (
        "cells.features",
        "cells.contours",
        "cells.expansion",
        "cells.classification",
        "layers.density",
        "layers.interface",
        "layers.support",
        "layers.arc",
        "layers.frame",
        "layers.band",
    ):
        importlib.import_module("granuscape." + name)


def test_unsupported_model_is_explicit():
    raw, masks = fixture()
    frame, _ = measure_cells(masks["soma"], masks["nuclei"], raw, 1)
    frame["analysis_ROI_fraction"] = 1.0
    classified, info = classify_cells(frame)
    assert not info["supported"]
    assert set(classified.morphotype) == {"unclassified_model_unsupported"}
    assert classified.mitral_probability.isna().all()


def test_pixel_area_not_conditional_on_object_ownership():
    soma = np.zeros((6, 6), np.uint32)
    soma[2:4, 2:4] = 1
    objects = np.zeros_like(soma)
    objects[1:4, 1:4] = 1
    table = object_table(
        objects, soma, np.ones_like(soma, bool), np.ones_like(soma), objects > 0, 1
    )
    assert table.soma_id.iloc[0] == 0
    assert np.count_nonzero((objects > 0) & (soma == 1)) == 4
    assert table.overlap_area_um2.iloc[0] == 9


def test_zero_puncta_retained_and_spatial_undefined():
    raw, masks = fixture()
    frame, _ = measure_cells(masks["soma"], masks["nuclei"], raw, 1)
    empty = np.zeros_like(masks["soma"])
    tables = {
        k: object_table(empty, masks["soma"], empty > 0, raw["af488"], empty > 0, 1)
        for k in ("g3bp", "eif3")
    }
    result = cell_metrics(
        frame, masks["soma"], masks["nuclei"], {"g3bp": empty, "eif3": empty}, tables, 1
    )
    assert len(result) == 5
    assert (result.g3bp_count == 0).all()
    assert (result.union_coverage == 0).all()
    assert result.g3bp_nearest_neighbor_um.isna().all()
    assert result.g3bp_angular_resultant.isna().all()


def test_refinement_protects_original_and_neighbors():
    raw, masks = fixture()
    profiles, _ = measure_cells(masks["soma"], masks["nuclei"], raw, 1)
    coverage = profiles.copy()
    coverage["analysis_ROI_fraction"] = 1.0
    labels, report = refine_cells(raw, masks, 1)
    assert np.array_equal(labels[masks["soma"] > 0], masks["soma"][masks["soma"] > 0])
    assert {r["soma_id"] for r in report} == set(profiles.soma_id)


def test_self_intersection_rejected():
    import pytest

    with pytest.raises(ValueError):
        validate_open_curve(np.array([[0, 0], [10, 10], [0, 10], [10, 0]], float))


def test_hdf5_tables_round_trip(tmp_path):
    frame = pd.DataFrame(
        {"id": [1, 2], "value": [0.0, np.nan], "status": ["zero", "undefined"]}
    )
    with h5py.File(tmp_path / "test.h5", "w") as h:
        write_table(h, "tables/cells", frame)
    with h5py.File(tmp_path / "test.h5", "r") as h:
        pd.testing.assert_frame_equal(frame, read_table(h, "tables/cells"))


def test_metadata_channel_identity_and_scale():
    xml = '<ImageDocument><Metadata><Information><Image><Dimensions><Channels><Channel Name="Rhoda"/><Channel Name="DAPI"/><Channel Name="AF488"/></Channels></Dimensions></Image></Information><Scaling><Items><Distance Id="X"><Value>0.000001</Value></Distance><Distance Id="Y"><Value>0.000001</Value></Distance></Items></Scaling></Metadata></ImageDocument>'
    mapping, scale = parse_metadata(
        xml, dict(dapi="DAPI", af488="AF488", rhoda="Rhoda")
    )
    assert mapping == dict(dapi=1, af488=2, rhoda=0)
    assert scale == 1


def test_integrated_arrays_analysis():
    raw, masks = fixture()
    curve = np.column_stack([np.full(256, 132.0), np.arange(256.0)])
    result = analyze_arrays(raw, masks, 1, boundary=curve)
    assert len(result["cells"]) == 5
    assert result["region"].any()
    assert result["image"]["review_required"]
    assert np.isclose(
        result["image"]["union_mcl_area_um2"],
        result["image"]["intersection_mcl_area_um2"],
    )
    assert not result["image"]["classification_supported"]


def test_expert_merge_and_frame_exclusion():
    import pytest
    from granuscape.cells.edits import apply_merges

    labels = np.zeros((20, 20), np.uint32)
    labels[5:10, 5:8] = 1
    labels[5:10, 8:11] = 2
    merged = apply_merges(labels, [dict(source_id=2, target_id=1)])
    assert np.count_nonzero(merged == 1) == 30
    assert np.count_nonzero(labels == 2) == 15
    labels[0:3, 0:3] = 3
    with pytest.raises(ValueError):
        apply_merges(labels, [dict(source_id=3, target_id=1)])


def test_automatic_layer_geometry():
    from granuscape.layers.detection import detect_layer

    raw, masks = fixture()
    y, x = np.indices(masks["nuclei"].shape)
    n = masks["nuclei"]
    sid = int(n.max()) + 1
    for cy in [0, 252]:
        for cx in range(145, 257, 12):
            n[(x - cx) ** 2 + (y - cy) ** 2 < 3**2] = sid
            sid += 1
    for cy in range(0, 256, 12):
        n[(x - 255) ** 2 + (y - cy) ** 2 < 3**2] = sid
        sid += 1
    curve, info = detect_layer(raw, masks, 1)
    assert curve.shape[1] == 2
    assert info["fit"]["fit_qc_passed"]
    assert np.isfinite(curve).all()


def test_batch_exports_and_isolates_failed_images(tmp_path, monkeypatch):
    from granuscape.pipeline import run as runner

    raw, masks = fixture()
    y, x = np.indices(masks["nuclei"].shape)
    n = masks["nuclei"]
    sid = int(n.max()) + 1
    for cy in [0, 252]:
        for cx in range(145, 257, 12):
            n[(x - cx) ** 2 + (y - cy) ** 2 < 3**2] = sid
            sid += 1
    for cy in range(0, 256, 12):
        n[(x - 255) ** 2 + (y - cy) ** 2 < 3**2] = sid
        sid += 1
    data = tmp_path / "data"
    data.mkdir()
    (data / "a.czi").write_bytes(b"fixture")
    (data / "b.czi").write_bytes(b"bad fixture")

    def imported(path, mapping):
        if path.stem == "b":
            raise ValueError("Intentional import failure")
        return raw, 1.0, "<fixture/>"

    monkeypatch.setattr(runner, "read_czi", imported)
    monkeypatch.setattr(
        runner,
        "infer",
        lambda raw, config, work, timer=None: (
            masks,
            {"backend": "synthetic test double"},
        ),
    )
    out = tmp_path / "results"
    config = dict(
        input=str(data),
        output=str(out),
        channels=dict(dapi="DAPI", af488="AF488", rhoda="Rhoda"),
        density_sigma_um=12,
    )
    assert runner.run_batch(config) == 1
    report = pd.read_csv(out / "run_summary.csv")
    assert list(report.status) == ["complete_requires_review", "failed"]
    assert report.loc[0, "elapsed_seconds"] >= report.loc[0, "analysis_seconds"] > 0
    assert (
        report.loc[0, "analysis_seconds"]
        >= report.loc[0, "analysis_soma_refinement_seconds"]
        > 0
    )
    assert report.loc[0, "qc_export_seconds"] > 0
    assert (out / "qc/a/cell_ids.png").is_file()
    assert not (out / "hdf5/a.partial.h5").exists()
    assert (out / "work/b/error.log").is_file()
    assert not (out / "work/a").exists()
    cells = pd.read_csv(out / "cells.csv")
    assert len(cells) == 5
    with h5py.File(out / "hdf5/a.h5", "r") as h:
        assert h.attrs["status"] == "complete_requires_review"
        assert "cells/labels_merged" in h
        assert len(read_table(h, "tables/cells")) == 5


def test_small_soma_expands_without_area_or_contrast_preselection():
    y, x = np.indices((64, 64))
    labels = np.zeros((64, 64), np.uint32)
    labels[30:33, 30:33] = 1
    nuclei = np.zeros_like(labels)
    nuclei[31, 31] = 1
    af = np.full(labels.shape, 10, np.uint16)
    af[(x - 31) ** 2 + (y - 31) ** 2 <= 2.5**2] = 100
    masks = dict(soma=labels, nuclei=nuclei)
    raw = dict(af488=af)
    refined, report = refine_cells(raw, masks, 1.0)
    entry = report[0]
    assert np.count_nonzero(labels == 1) == 9  # Below the former 45 um² cutoff.
    assert entry["signal_sample_pixels"] == 8  # Below the former 20-pixel cutoff.
    assert entry["status"] == "expanded_requires_review"
    assert entry["applied"]
    assert np.count_nonzero(refined == 1) > 9
    assert np.array_equal(refined[labels > 0], labels[labels > 0])


def test_every_seed_has_an_explicit_refinement_decision():
    raw, masks = fixture()
    # Zero apparent contrast does not pre-exclude a cell: the signal check records why it cannot grow.
    raw["af488"][:] = 20
    _, report = refine_cells(raw, masks, 1.0)
    assert len(report) == 5
    assert all(r["status"] == "insufficient_signal_review" for r in report)
    assert all(not r["applied"] for r in report)
