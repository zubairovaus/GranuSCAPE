"""Protect source data and the agreed QC export contract."""

import hashlib
import json
import h5py
import numpy as np
from PIL import Image
import pytest
from test_scientific_contract import fixture
from granuscape.cells.refinement import measure_cells
from granuscape.io.hdf5 import write_table
from granuscape.qc.export import export_qc, load_display_config
from granuscape.qc.plots import save_qc, composite
from granuscape.qc.settings import resolve_options


def data():
    raw, masks = fixture()
    cells, _ = measure_cells(masks["soma"], masks["nuclei"], raw, 1.0)
    cells["analysis_ROI_fraction"] = [1, 1, 1, 0, 0]
    curve = np.array([[132.0, 0.0], [132.0, 255.0]])
    roi = np.zeros(masks["soma"].shape, bool)
    roi[:, 110:150] = True
    return raw, masks, cells, curve, roi


def test_default_views_keep_ids_and_do_not_modify_arrays(tmp_path):
    raw, masks, cells, curve, roi = data()
    originals = {k: v.copy() for k, v in dict(raw, **masks).items()}
    outputs = save_qc(
        raw, masks["soma"], masks["nuclei"], masks, curve, roi, cells, 1.0, tmp_path
    )
    assert set(outputs) == {
        "A_overlay.png",
        "B_masks.png",
        "C_cells.png",
        "D_puncta.png",
        "cell_ids.png",
    }
    for key, value in dict(raw, **masks).items():
        assert np.array_equal(value, originals[key])
    for filename in outputs:
        with Image.open(tmp_path / filename) as image:
            image.verify()
    with Image.open(tmp_path / "A_overlay.png") as a, Image.open(
        tmp_path / "cell_ids.png"
    ) as numbered:
        assert a.size == numbered.size
        assert a.width == masks["soma"].shape[1]
        assert not np.array_equal(np.array(a), np.array(numbered))
    settings = json.loads((tmp_path / "display_settings.json").read_text())
    assert settings["options"]["mcl_color"] == "#C86F80"
    assert settings["options"]["mcl_line_width_px"] == 3


def test_selected_views_and_optional_id_only(tmp_path):
    raw, masks, cells, curve, roi = data()
    output = save_qc(
        raw,
        masks["soma"],
        masks["nuclei"],
        masks,
        curve,
        roi,
        cells,
        1.0,
        tmp_path / "masks",
        dict(views=["B", "D"], include_cell_ids=False),
    )
    assert set(output) == {"B_masks.png", "D_puncta.png"}
    output = save_qc(
        raw,
        masks["soma"],
        masks["nuclei"],
        masks,
        curve,
        roi,
        cells,
        1.0,
        tmp_path / "ids",
        dict(views=[], include_cell_ids=True),
    )
    assert output == ["cell_ids.png"]


def test_invalid_display_configuration_fails():
    for options in (
        {"views": ["E"]},
        {"views": ["A", "A"]},
        {"include_cell_ids": "yes"},
        {"mcl_color": "red"},
        {"mcl_line_width_px": 0},
        {"clahe_kernel_px": [0, 128]},
        {"normalization": "unknown"},
        {"clahe_clip_limit": float("nan")},
        {"typo": 1},
    ):
        with pytest.raises(ValueError):
            resolve_options(options)


def test_qc_configuration_does_not_require_backends(tmp_path):
    config = tmp_path / "qc.yaml"
    config.write_text("qc:\n  views: [C]\n  include_cell_ids: false\n")
    result = load_display_config(config)
    assert result["views"] == ["C"]
    assert result["normalization"] == "clahe"


def test_constant_channels_stay_black():
    raw = {k: np.full((32, 32), 12, np.uint16) for k in ("rhoda", "af488", "dapi")}
    assert not composite(raw).any()


def test_readonly_regeneration_and_existing_output_protection(tmp_path):
    raw, masks, cells, curve, roi = data()
    source = tmp_path / "image.h5"
    with h5py.File(source, "w") as h:
        h.attrs.update(schema="granuscape", status="complete_requires_review")
        for key, value in raw.items():
            h.create_dataset("raw/" + key, data=value)
        h.create_dataset("cells/labels_refined", data=masks["soma"])
        h.create_dataset("segmentation/nuclei", data=masks["nuclei"])
        for key in ("g3bp", "eif3"):
            h.create_dataset("segmentation/" + key, data=masks[key])
        h.create_dataset("mcl/boundary_xy_px", data=curve)
        h.create_dataset("mcl/region", data=roi)
        h.create_group("acquisition").attrs["pixel_size_um"] = 1.0
        write_table(h, "tables/cells", cells)
    before = hashlib.sha256(source.read_bytes()).hexdigest()
    assert (
        export_qc(source, tmp_path / "qc", dict(views=["B"], include_cell_ids=True))
        == 0
    )
    assert hashlib.sha256(source.read_bytes()).hexdigest() == before
    assert (tmp_path / "qc/image/B_masks.png").is_file()
    assert (tmp_path / "qc/image/cell_ids.png").is_file()
    with pytest.raises(FileExistsError):
        export_qc(source, tmp_path / "qc")
