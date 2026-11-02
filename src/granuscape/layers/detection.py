"""Anatomical layer interface and adaptive nuclear-density band."""

import numpy as np
from .density import nucleus_density
from .interface import choose_sparse_interface, validate_open_curve
from .support import soma_features, fit_density_af488
from .frame import extend_spline
from .band import compute_band


def detect_layer(raw, masks, pixel_size, sigma_um=12):
    density = nucleus_density(masks["nuclei"], sigma_um, pixel_size)
    interface, interface_info = choose_sparse_interface(
        density, pixel_size, gap_closing_um=sigma_um
    )
    features = soma_features(masks["soma"], masks["nuclei"], raw["af488"], pixel_size)
    model, fit_info, _ = fit_density_af488(interface * pixel_size, features, pixel_size)
    _, curve, curvature, _ = extend_spline(model, pixel_size, masks["soma"].shape)
    validate_open_curve(curve)
    if np.max(abs(curvature)) > 1 / 20:
        raise ValueError(
            "MCL curvature exceeds the supported anatomical range; review required"
        )
    return curve, {"interface": interface_info, "fit": fit_info}


def adaptive_band(raw, labels, nuclei, pixel_size, curve, profiles):
    return compute_band(
        labels,
        nuclei,
        raw["af488"],
        raw["dapi"],
        pixel_size,
        curve,
        profiles.to_dict("records"),
    )
