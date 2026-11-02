"""Analyze one set of native channels and segmentation arrays."""

from granuscape.benchmark.timing import Timer
import numpy as np
from granuscape.cells.refinement import measure_cells, refine_cells, classify_cells
from granuscape.layers.detection import detect_layer, adaptive_band
from granuscape.granules.objects import object_table
from granuscape.measurements.profiles import cell_metrics


def validate_arrays(raw, masks, pixel_size):
    shape = raw["dapi"].shape
    if len(shape) != 2 or pixel_size <= 0 or (not np.isfinite(pixel_size)):
        raise ValueError("Expected calibrated two-dimensional native planes")
    for name in ("dapi", "af488", "rhoda"):
        if raw[name].shape != shape or not np.isfinite(raw[name]).all():
            raise ValueError("Invalid raw channel: " + name)
    for name in ("soma", "nuclei", "g3bp", "eif3"):
        array = masks[name]
        if array.shape != shape or array.dtype.kind not in "iu" or array.min() < 0:
            raise ValueError("Invalid segmentation: " + name)
    if not masks["soma"].any() or not masks["nuclei"].any():
        raise ValueError("No soma or nuclear segmentation; image requires review")


def analyze_arrays(
    raw,
    masks,
    pixel_size,
    sigma_um=12,
    boundary=None,
    edits=None,
    timer=None,
    spatial=None,
    image_id="",
):
    timer = Timer() if timer is None else timer
    validate_arrays(raw, masks, pixel_size)
    from granuscape.cells.edits import apply_merges

    masks = dict(masks)
    masks["soma"] = apply_merges(masks["soma"], edits or [])
    merged = masks["soma"].copy()
    if boundary is None:
        with timer.span("analysis_layer_detection"):
            boundary, layer_info = detect_layer(raw, masks, pixel_size, sigma_um)
    else:
        from granuscape.layers.interface import validate_open_curve

        validate_open_curve(boundary)
        if (
            np.any(boundary < 0)
            or np.any(boundary[:, 0] > masks["soma"].shape[1] - 1)
            or np.any(boundary[:, 1] > masks["soma"].shape[0] - 1)
        ):
            raise ValueError("Explicit MCL boundary leaves the native image grid")
        layer_info = {"source": "explicit_boundary", "review_required": True}
    import pandas as pd

    with timer.span("analysis_soma_refinement"):
        refined, expansion = refine_cells(raw, masks, pixel_size)
    with timer.span("analysis_soma_measurements"):
        profiles, ranks = measure_cells(refined, masks["nuclei"], raw, pixel_size)
    with timer.span("analysis_adaptive_ROI"):
        roi, core, coverage, band, band_info = adaptive_band(
            raw, refined, masks["nuclei"], pixel_size, boundary, profiles
        )
    coverage = pd.DataFrame(coverage)[
        ["soma_id", "analysis_ROI_fraction", "density_core_fraction", "soma_support"]
    ]
    profiles = profiles.merge(coverage, on="soma_id", validate="one_to_one")
    reports = {r["soma_id"]: r for r in expansion}
    profiles["expansion_status"] = [reports[sid]["status"] for sid in profiles.soma_id]
    profiles["expansion_applied"] = [
        reports[sid].get("applied", False) for sid in profiles.soma_id
    ]
    profiles["expansion_review_flags"] = [
        ";".join(reports[sid].get("review_flags", [])) for sid in profiles.soma_id
    ]
    profiles["expansion_error"] = [
        reports[sid].get("error", "") for sid in profiles.soma_id
    ]
    profiles["expansion_signal_pixels"] = [
        reports[sid].get("signal_sample_pixels", np.nan) for sid in profiles.soma_id
    ]
    original_areas = np.bincount(merged.ravel()) * pixel_size**2
    profiles["soma_seed_area_um2"] = [original_areas[sid] for sid in profiles.soma_id]
    unstable = {
        r["soma_id"]
        for r in expansion
        if r.get("review_flags") or ("review" in r["status"] and (not r.get("applied")))
    }
    profiles.loc[profiles.soma_id.isin(unstable), "measurement_valid"] = False
    layer_fit_passed = bool(layer_info.get("fit", {}).get("fit_qc_passed", True))
    if not layer_fit_passed:
        profiles["measurement_valid"] = False
    with timer.span("analysis_classification"):
        profiles, classification = classify_cells(profiles)
    tables = {}
    for marker, channel, other in (
        ("g3bp", "rhoda", "eif3"),
        ("eif3", "af488", "g3bp"),
    ):
        with timer.span("analysis_puncta_tables"):
            tables[marker] = object_table(
                masks[marker], refined, roi, raw[channel], masks[other] > 0, pixel_size
            )
    with timer.span("analysis_spatial_profiles"):
        profiles = cell_metrics(
            profiles, refined, masks["nuclei"], masks, tables, pixel_size
        )
    from granuscape.spatial import spatial_analysis

    with timer.span("analysis_projection_spatial"):
        spatial_result = spatial_analysis(
            profiles,
            refined,
            masks["nuclei"],
            masks,
            tables,
            pixel_size,
            spatial,
            image_id,
        )
    tables = spatial_result["objects"]
    image = {
        "mcl_area_um2": float(roi.sum() * pixel_size**2),
        "nuclei_intersecting_mcl": int(
            len(np.unique(masks["nuclei"][roi]))
            - int(np.any(masks["nuclei"][roi] == 0))
        ),
        "cells_intersecting_mcl": int((profiles.analysis_ROI_fraction >= 0.1).sum()),
        "mitral_like_candidates": int(
            (
                (profiles.morphotype == "mitral_like")
                & (profiles.analysis_ROI_fraction >= 0.1)
            ).sum()
        ),
        "classification_supported": bool(classification["supported"]),
        "review_required": True,
        "mcl_fit_qc_passed": layer_fit_passed,
    }
    a, b = (masks["g3bp"] > 0, masks["eif3"] > 0)
    for marker, positive in (
        ("g3bp", a),
        ("eif3", b),
        ("intersection", a & b),
        ("union", a | b),
    ):
        image[marker + "_mcl_area_um2"] = float((positive & roi).sum() * pixel_size**2)
        image[marker + "_mcl_coverage"] = (
            float(positive[roi].mean()) if roi.any() else np.nan
        )
        if marker in tables:
            image[marker + "_mcl_object_count"] = int(
                (tables[marker].in_mcl_fraction > 0).sum()
            )
    return dict(
        spatial=spatial_result,
        labels_merged=merged,
        labels_refined=refined,
        boundary=boundary,
        region=roi,
        core=core,
        band=band,
        cells=profiles,
        objects=tables,
        image=image,
        diagnostics={
            "layer": layer_info,
            "band": band_info,
            "classification": classification,
            "expansion": expansion,
            "nucleus_associations": ranks,
        },
    )
