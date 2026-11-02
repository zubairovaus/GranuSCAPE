import numpy as np
from scipy.optimize import brentq
from .arc import curvature_sign_changes


def extend_spline(model, px, shape, n_original=1600):
    height, width = shape

    def margin(t):
        x, y = model(t) / px
        return min(x, width - 1 - x, y, height - 1 - y)

    roots = []
    if margin(0.0) <= 0 or margin(1.0) <= 0:
        grid = np.linspace(0.0, 1.0, max(4097, n_original * 2 + 1))
        xy = model(grid) / px
        margins = np.min(
            np.c_[xy[:, 0], width - 1 - xy[:, 0], xy[:, 1], height - 1 - xy[:, 1]],
            axis=1,
        )
        inside = margins > 1e-08
        transitions = np.diff(np.r_[False, inside, False].astype(int))
        begins = np.flatnonzero(transitions == 1)
        ends = np.flatnonzero(transitions == -1) - 1
        if len(begins) != 1:
            raise ValueError(
                "Arc has no unique connected in-frame component; manual review required"
            )
        a, b = (int(begins[0]), int(ends[0]))
        for endpoint, direction, index in [(0.0, -1.0, a), (1.0, 1.0, b)]:
            if abs(margin(endpoint)) <= 1e-08:
                root = endpoint
            elif margin(endpoint) < 0:
                lo, hi = (
                    (grid[index - 1], grid[index])
                    if direction < 0
                    else (grid[index], grid[index + 1])
                )
                root = brentq(margin, lo, hi, xtol=1e-14)
            else:
                previous = endpoint
                for i in range(1, 201):
                    trial = endpoint + direction * i * 0.01
                    if margin(trial) <= 0:
                        root = brentq(
                            margin,
                            min(previous, trial),
                            max(previous, trial),
                            xtol=1e-14,
                        )
                        break
                    previous = trial
                else:
                    raise ValueError("No nearby frame crossing: manual review required")
            roots.append(root)
        start, end = roots
        if end <= start:
            raise ValueError("Invalid in-frame parameter interval")
        t = np.linspace(start, end, n_original + 398)
        curve = model(t) / px
        if not np.isfinite(curve).all():
            raise ValueError("Nonfinite frame-intersected arc")
        if (
            np.any(curve < -1e-06)
            or np.any(curve[:, 0] > width - 1 + 1e-06)
            or np.any(curve[:, 1] > height - 1 + 1e-06)
        ):
            raise ValueError(
                "Frame-intersected curve leaves frame; manual review required"
            )
        curve[:, 0] = np.clip(curve[:, 0], 0, width - 1)
        curve[:, 1] = np.clip(curve[:, 1], 0, height - 1)
        d1 = model.derivative(1)(t)
        d2 = model.derivative(2)(t)
        if np.any(np.linalg.norm(d1, axis=1) <= 1e-12):
            raise ValueError("Degenerate arc derivative")
        k = (d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]) / np.linalg.norm(
            d1, axis=1
        ) ** 3
        if not curvature_sign_changes(k) == 0:
            raise ValueError("Frame-intersected curve changed bend direction")
        return (t, curve, k, roots)
    for endpoint, direction in [(0.0, -1.0), (1.0, 1.0)]:
        previous = endpoint
        if margin(endpoint) <= 0:
            raise ValueError("Existing endpoint must be inside frame")
        for i in range(1, 201):
            trial = endpoint + direction * i * 0.01
            if margin(trial) <= 0:
                roots.append(
                    brentq(
                        margin, min(previous, trial), max(previous, trial), xtol=1e-14
                    )
                )
                break
            previous = trial
        else:
            raise ValueError("No nearby frame crossing: manual review required")
    start, end = roots
    t = np.r_[
        np.linspace(start, 0, 200)[:-1],
        np.linspace(0, 1, n_original),
        np.linspace(1, end, 200)[1:],
    ]
    curve = model(t) / px
    if not (np.all(curve[:, 0] >= -1e-06) and np.all(curve[:, 0] <= width - 1 + 1e-06)):
        raise ValueError("Numerical invariant violated")
    if not (
        np.all(curve[:, 1] >= -1e-06) and np.all(curve[:, 1] <= height - 1 + 1e-06)
    ):
        raise ValueError("Numerical invariant violated")
    curve[:, 0] = np.clip(curve[:, 0], 0, width - 1)
    curve[:, 1] = np.clip(curve[:, 1], 0, height - 1)
    d1 = model.derivative(1)(t)
    d2 = model.derivative(2)(t)
    k = (d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]) / np.linalg.norm(d1, axis=1) ** 3
    if not curvature_sign_changes(k) == 0:
        raise ValueError("Extended curve changed bend direction")
    return (t, curve, k, roots)
