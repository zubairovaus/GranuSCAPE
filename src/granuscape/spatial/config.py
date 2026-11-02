"""Explicit spatial model configuration; Monte Carlo is disabled by default."""

import numpy as np

DEFAULTS = dict(
    perinuclear_radius_um=2.0,
    monte_carlo=False,
    simulations=1000,
    seed=17,
    modes=["full_mask", "centroid_only"],
    max_candidates=200000,
    diagnostics=True,
    diagnostic_cells=3,
)


def options(value=None):
    value = {} if value is None else dict(value)
    if set(value) - set(DEFAULTS):
        raise ValueError(
            "Unknown spatial options: " + str(sorted(set(value) - set(DEFAULTS)))
        )
    result = {**DEFAULTS, **value}
    for k in ("monte_carlo", "diagnostics"):
        if not isinstance(result[k], bool):
            raise ValueError(k + " must be boolean")
    for k in ("simulations", "max_candidates", "diagnostic_cells", "seed"):
        if (
            isinstance(result[k], bool)
            or not isinstance(result[k], int)
            or result[k] < (1 if k in ("simulations", "max_candidates") else 0)
        ):
            raise ValueError(k + " must be a nonnegative integer (positive for counts)")
    if (
        not np.isfinite(result["perinuclear_radius_um"])
        or result["perinuclear_radius_um"] <= 0
    ):
        raise ValueError("Invalid perinuclear radius")
    if (
        not isinstance(result["modes"], list)
        or not result["modes"]
        or set(result["modes"]) - {"full_mask", "centroid_only"}
        or len(set(result["modes"])) != len(result["modes"])
    ):
        raise ValueError("Invalid Monte Carlo modes")
    return result
