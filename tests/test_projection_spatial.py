import numpy as np, pandas as pd, h5py, pytest
from granuscape.spatial.geometry import regions, angular_domains, field_metrics
from granuscape.spatial.nulls import valid_anchors, simulate, summarize, resultant
from granuscape.spatial.config import options
from granuscape.spatial.analysis import spatial_analysis
from granuscape.granules.objects import object_table
from granuscape.spatial.export import reanalyze
from granuscape.io.hdf5 import write_table, read_table


def fixture():
    s = np.zeros((21, 21), np.int32)
    s[2:19, 2:19] = 1
    n = np.zeros_like(s)
    n[8:13, 8:13] = 7
    n[3:5, 3:5] = 8
    g = np.zeros_like(s)
    g[10, 10] = 1
    g[10, 12:15] = 2
    g[4, 4] = 3
    g[17, 10] = 4
    e = np.zeros_like(s)
    e[6, 10] = 1
    c = pd.DataFrame([dict(soma_id=1, primary_nucleus_id=7)])
    m = dict(g3bp=g, eif3=e)
    o = {
        k: object_table(v, s, s > 0, np.zeros(s.shape), np.zeros(s.shape, bool), 1.0)
        for k, v in m.items()
    }
    return c, s, n, m, o


def test_partition_primary_priority_and_foreign():
    s = np.ones((7, 7), bool)
    p = np.zeros_like(s)
    p[3, 3] = 1
    a = p.copy()
    a[1, 1] = 1
    r = regions(s, p, a, 1.0)
    assert [r[k].sum() for k in ["nuclear", "foreign", "perinuclear", "distal"]] == [
        1,
        1,
        12,
        35,
    ]
    assert r["associated"].sum() == 13
    assert np.array_equal(
        sum(r[k].astype(int) for k in ["nuclear", "foreign", "perinuclear", "distal"]),
        s.astype(int),
    )
    s[0, :] = False
    p[:2, 2:4] = True
    r = regions(s, p, a, 0.5)
    assert np.array_equal(r["nuclear"], s & p)


def test_objects_and_no_double_count():
    c, s, n, m, o = fixture()
    r = spatial_analysis(c, s, n, m, o, 1.0, {"diagnostics": False})
    a = r["object_regions"].query('marker=="g3bp"').set_index("punctum_id")
    row = r["cells"].query('marker=="g3bp"').iloc[0]
    assert (
        len(r["objects"]["g3bp"]) == 4
        and a.loc[1, "spatial_class"] == "nuclear_overlapping"
    )
    assert a.loc[2, "nuclear_fraction_of_whole_object"] == pytest.approx(1 / 3)
    assert a.loc[2, "perinuclear_fraction_of_whole_object"] == pytest.approx(2 / 3)
    assert a.loc[3, "spatial_class"] == "ambiguous_foreign_nucleus_associated"
    assert (
        row.associated_owned_object_count == 2
        and row.nuclear_owned_object_count + row.perinuclear_owned_object_count == 3
    )
    assert row.foreign_positive_area_um2 == 1 and row.nuclear_positive_area_um2 == 2
    assert a.loc[1, "centroid_signed_distance_um"] < 0


def test_foreign_primary_spanning_and_soma_border():
    c, s, n, m, o = fixture()
    for coords in [(slice(4, 11), slice(4, 11)), (slice(0, 4), 10)]:
        m["g3bp"][:] = 0
        m["g3bp"][coords] = 9
        o["g3bp"] = object_table(
            m["g3bp"], s, s > 0, np.zeros(s.shape), np.zeros(s.shape, bool), 1.0
        )
        a = (
            spatial_analysis(c, s, n, m, o, 1.0, {"diagnostics": False})[
                "object_regions"
            ]
            .query('marker=="g3bp"')
            .iloc[0]
        )
        if isinstance(coords[1], slice):
            assert (
                a.nuclear_overlap_area_um2 > 0
                and a.foreign_overlap_area_um2 > 0
                and a.assignment_ambiguous
            )
        else:
            assert a.object_area_um2 == 4 and a.soma_intersection_area_um2 == 2


def test_empty_and_zero_policy():
    s = np.ones((5, 5), bool)
    p = s.copy()
    r = regions(s, p, p, 1.0)
    sec, _ = angular_domains(s.shape, p)
    v = field_metrics(s, r, sec)
    assert (
        np.isnan(v["perinuclear_coverage"])
        and np.isnan(v["distal_coverage"])
        and v["nuclear_coverage"] == 1
    )
    p = np.zeros_like(s)
    p[2, 2] = 1
    r = regions(s, p, p, 0.25)
    sec, _ = angular_domains(s.shape, p)
    v = field_metrics(np.zeros_like(s), r, sec)
    assert v["perinuclear_coverage"] == 0 and np.isnan(v["perinuclear_enrichment"])
    assert np.isnan(resultant(np.array([[1.0, 1.0]]), (2, 2)))


def test_shape_placement_and_failure():
    a = np.ones((8, 9), bool)
    a[2:5, 2:5] = False
    c = np.array([[0, 0], [1, 0], [1, 1]])
    anchors, error = valid_anchors(a, c, 1000)
    assert not error
    for anchor in anchors:
        assert a[tuple((c + anchor).T)].all()
    assert not len(valid_anchors(np.zeros((3, 3), bool), c, 100)[0])
    assert (
        valid_anchors(np.ones((100, 100), bool), c, 2)[1] == "placement_candidate_limit"
    )


def test_simulation_reproducibility_and_domains():
    s = np.ones((11, 11), bool)
    p = np.zeros_like(s)
    p[5, 5] = 1
    r = regions(s, p, p, 1)
    sec, center = angular_domains(s.shape, p)
    pos = np.zeros_like(s)
    pos[2, 2:4] = 1
    obj = [dict(coords=np.argwhere(pos))]
    opt = options({"simulations": 30})
    a, v, first = simulate(obj, pos, r, sec, center, opt, 2, "strict", "full_mask")
    b, v2, _ = simulate(obj, pos, r, sec, center, opt, 2, "strict", "full_mask")
    assert all(np.array_equal(v[k], v2[k], equal_nan=True) for k in v)
    assert first.sum() == 2 and not first[p].any()
    assert a[0]["successful"] == 30 and 0 < a[0]["empirical_p"] <= 1
    a, _, _ = simulate(
        [dict(coords=np.argwhere(p))], p, r, sec, center, opt, 2, "strict", "full_mask"
    )
    assert "crosses_allowed" in a[0]["model_qc"]
    a, _, first = simulate(
        [dict(coords=np.argwhere(p))],
        p,
        r,
        sec,
        center,
        opt,
        2,
        "extended",
        "full_mask",
    )
    assert first.sum() == 1


def test_centroid_asymmetry_fragmentation():
    s = np.zeros((15, 15), bool)
    s[1:6, 1:5] = 1
    s[9:14, 8:14] = 1
    p = np.zeros_like(s)
    p[3, 3] = 1
    r = regions(s, p, p, 1)
    sec, center = angular_domains(s.shape, p)
    coords = np.array([[2, 2], [12, 12]])
    pos = np.zeros_like(s)
    pos[tuple(coords.T)] = 1
    a, v, _ = simulate(
        [dict(coords=c[None, :]) for c in coords],
        pos,
        r,
        sec,
        center,
        options({"simulations": 20}),
        7,
        "strict",
        "centroid_only",
    )
    assert a[0]["successful"] == 20


def test_statistics():
    a = summarize(2.0, [0, 0, 0], 3)
    assert (
        a["empirical_p"] == 0.5
        and a["qc"] == "zero_null_variance"
        and np.isnan(a["z_score"])
    )
    assert summarize(np.nan, [1, 2], 2)["qc"] == "observed_undefined"


def test_readonly_export(tmp_path):
    c, s, n, m, o = fixture()
    path = tmp_path / "test.h5"
    with h5py.File(path, "w") as h:
        h.attrs["schema"] = "granuscape"
        h.create_dataset("cells/labels_refined", data=s)
        h.create_dataset("segmentation/nuclei", data=n)
        h.require_group("acquisition").attrs["pixel_size_um"] = 1.0
        write_table(h, "tables/cells", c)
        for marker in m:
            h.create_dataset("segmentation/" + marker, data=m[marker])
            write_table(h, "tables/puncta_" + marker, o[marker])
    before = path.read_bytes()
    assert (
        reanalyze(
            path,
            tmp_path / "out",
            {"simulations": 5, "monte_carlo": True, "diagnostics": False},
        )
        == 0
    )
    assert before == path.read_bytes()
    with h5py.File(tmp_path / "out/test_spatial.h5") as h:
        assert (
            len(read_table(h, "tables/spatial_cells")) == 4 and "null" in h["spatial"]
        )


def test_invalid_options():
    for value in [
        {"simulations": 0},
        {"seed": True},
        {"monte_carlo": "yes"},
        {"perinuclear_radius_um": -1},
        {"modes": ["pixel"]},
        {"unused": 1},
    ]:
        with pytest.raises(ValueError):
            options(value)


def test_primary_missing_and_merged_nucleus_qc():
    c, s, n, m, o = fixture()
    c.primary_nucleus_id = 0
    r = spatial_analysis(c, s, n, m, o, 1.0, {"diagnostics": False})
    assert (
        len(r["cells"]) == 4
        and r["cells"].spatial_qc.eq("primary_nucleus_missing").all()
    )
    c.primary_nucleus_id = 7
    n[3:5, 10:12] = 7
    r = spatial_analysis(c, s, n, m, o, 1.0, {"diagnostics": False})
    assert r["cells"].spatial_qc.str.contains("disconnected_primary_nucleus").all()


def test_radial_zones_reconstruct_strict_area_and_coverage():
    c, s, n, m, o = fixture()
    r = spatial_analysis(c, s, n, m, o, 1.0, {"diagnostics": False})
    rad = r["radial_profiles"].query('marker=="g3bp"')
    row = r["cells"].query('marker=="g3bp"').iloc[0]
    assert rad.zone_area_um2.sum() == row.perinuclear_area_um2 + row.distal_area_um2
    peri = rad[rad.distance_high_um <= 2]
    assert peri.positive_area_um2.sum() / peri.zone_area_um2.sum() == pytest.approx(
        row.perinuclear_coverage
    )


def test_no_objects_yields_explicit_null_qc():
    s = np.ones((7, 7), bool)
    p = np.zeros_like(s)
    p[3, 3] = 1
    r = regions(s, p, p, 1)
    sec, center = angular_domains(s.shape, p)
    stats, arr, first = simulate(
        [],
        np.zeros_like(s),
        r,
        sec,
        center,
        options({"simulations": 10}),
        17,
        "strict",
        "full_mask",
    )
    assert first is None and all(
        x["successful"] == 0 and "no_objects" in x["model_qc"] for x in stats
    )


def test_pipeline_adds_spatial_without_changing_legacy_metrics():
    from test_scientific_contract import fixture as pipeline_fixture
    from granuscape.pipeline.analysis import analyze_arrays

    raw, masks = pipeline_fixture()
    boundary = np.array([[130.0, 0.0], [130.0, 128.0], [130.0, 255.0]])
    a = analyze_arrays(
        raw, masks, 1.0, boundary=boundary, spatial={"diagnostics": False}
    )
    b = analyze_arrays(
        raw,
        masks,
        1.0,
        boundary=boundary,
        spatial={"monte_carlo": True, "simulations": 5, "diagnostics": False},
    )
    pd.testing.assert_frame_equal(a["cells"], b["cells"])
    assert len(a["spatial"]["cells"]) == 20
    assert len(b["spatial"]["null_summary"]) > 0
    assert np.array_equal(a["labels_refined"], b["labels_refined"])
