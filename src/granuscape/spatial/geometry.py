"""Disjoint two-dimensional analysis regions with primary-nucleus priority."""

import numpy as np
from scipy import ndimage as ndi


def regions(soma, primary, all_nuclei, pixel_size, radius=2.0):
    soma, primary, all_nuclei = map(
        lambda a: np.asarray(a, bool), (soma, primary, all_nuclei)
    )
    if not (soma.shape == primary.shape == all_nuclei.shape):
        raise ValueError("Spatial masks must share the native grid")
    if (
        not np.isfinite(pixel_size)
        or pixel_size <= 0
        or not np.isfinite(radius)
        or radius <= 0
    ):
        raise ValueError("Positive finite calibration and radius required")
    if not primary.any():
        raise ValueError("Primary nucleus is absent")
    all_nuclei = all_nuclei | primary
    distance = ndi.distance_transform_edt(~primary, sampling=pixel_size)
    nuclear = soma & primary
    foreign = soma & all_nuclei & ~primary
    allowed = soma & ~all_nuclei
    peri = allowed & (distance <= radius)
    distal = allowed & (distance > radius)
    assert np.array_equal(
        nuclear.astype(int) + foreign + peri + distal, soma.astype(int)
    )
    return dict(
        nuclear=nuclear,
        foreign=foreign,
        perinuclear=peri,
        distal=distal,
        associated=nuclear | peri,
        strict=allowed,
        extended=allowed | nuclear,
        distance=distance,
        signed_distance=distance
        - ndi.distance_transform_edt(primary, sampling=pixel_size),
    )


def angular_domains(shape, primary):
    cy, cx = ndi.center_of_mass(primary)
    yy, xx = np.indices(shape)
    angle = np.arctan2(yy - cy, xx - cx)
    sectors = np.minimum(7, np.floor((angle + np.pi) / (np.pi / 4)).astype(int))
    return sectors, (cy, cx)


def field_metrics(positive, region, sectors):
    result = {}
    for name in ("nuclear", "foreign", "perinuclear", "distal", "associated"):
        zone = region[name]
        result[name + "_coverage"] = (
            float(positive[zone].mean()) if zone.any() else np.nan
        )
    distal = result["distal_coverage"]
    for prefix, name in [
        ("perinuclear", "perinuclear"),
        ("nucleus_associated", "associated"),
    ]:
        value = result[name + "_coverage"]
        result[prefix + "_enrichment"] = (
            value / distal if np.isfinite(distal) and distal > 0 else np.nan
        )
        result[prefix + "_coverage_difference"] = value - distal
    for variant in ("strict", "extended"):
        values = np.array(
            [
                (
                    positive[region[variant] & (sectors == i)].mean()
                    if (region[variant] & (sectors == i)).any()
                    else np.nan
                )
                for i in range(8)
            ]
        )
        result[variant + "_angular_coverage_cv"] = (
            float(values.std() / values.mean())
            if np.isfinite(values).all() and values.mean() > 0
            else np.nan
        )
    return result
