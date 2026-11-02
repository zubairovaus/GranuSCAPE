import numpy as np
import pandas as pd
from scipy.ndimage import find_objects, distance_transform_edt
from scipy.spatial import cKDTree
from scipy.interpolate import BSpline
from scipy.optimize import minimize, LinearConstraint

SETTINGS = dict(
    soma_area_min_um2=45.0,
    ring_inner_um=3.0,
    ring_outer_um=12.0,
    min_cytoplasm_pixels=20,
    min_background_pixels=30,
    corridor_um=20.0,
    af488_weight_per_bin=0.35,
    bins=40,
    smoothing_penalty=0.01,
    max_shift_from_density_um=10.0,
)


def soma_features(soma, nuc, af488, px):
    if soma.shape != nuc.shape or soma.shape != af488.shape or soma.ndim != 2:
        raise ValueError("Expected matching 2D soma, nuclei and native AF488")
    if not np.isfinite(af488).all() or px <= 0:
        raise ValueError("Invalid raw intensity/scale")
    pad = int(np.ceil(SETTINGS["ring_outer_um"] / px)) + 1
    rows = []
    for label, sl in enumerate(find_objects(soma), 1):
        if sl is None:
            continue
        y0 = max(0, sl[0].start - pad)
        y1 = min(soma.shape[0], sl[0].stop + pad)
        x0 = max(0, sl[1].start - pad)
        x1 = min(soma.shape[1], sl[1].stop + pad)
        labels = soma[y0:y1, x0:x1]
        n = nuc[y0:y1, x0:x1]
        a = af488[y0:y1, x0:x1]
        mask = labels == label
        ys, xs = np.nonzero(mask)
        area = float(mask.sum() * px * px)
        row = dict(
            soma_id=label,
            x_px=float(xs.mean() + x0),
            y_px=float(ys.mean() + y0),
            area_um2=area,
            af488_cytoplasm_median=np.nan,
            af488_background_median=np.nan,
            local_contrast_z=0.0,
            brightness_weight=0.0,
            feature_valid=False,
        )
        if area >= SETTINGS["soma_area_min_um2"]:
            dist = distance_transform_edt(~mask) * px
            ring = (
                (dist >= SETTINGS["ring_inner_um"])
                & (dist <= SETTINGS["ring_outer_um"])
                & (labels == 0)
                & (n == 0)
            )
            cyto = mask & (n == 0)
            if (
                cyto.sum() >= SETTINGS["min_cytoplasm_pixels"]
                and ring.sum() >= SETTINGS["min_background_pixels"]
            ):
                signal = float(np.median(a[cyto]))
                bg = float(np.median(a[ring]))
                mad = float(np.median(np.abs(a[ring].astype(float) - bg)))
                scale = max(1.4826 * mad, abs(bg) * 0.01, np.finfo(float).eps)
                z = max(0.0, (signal - bg) / scale)
                row.update(
                    af488_cytoplasm_median=signal,
                    af488_background_median=bg,
                    local_contrast_z=z,
                    brightness_weight=z / (1.0 + z),
                    feature_valid=True,
                )
        rows.append(row)
    return pd.DataFrame(rows)


def fit_density_af488(raw_um, features, px):
    raw = np.asarray(raw_um, float)
    origin = raw.mean(0)
    _, _, vh = np.linalg.svd(raw - origin, full_matrices=False)
    along, across = vh
    s = (raw - origin) @ along
    if s[-1] < s[0]:
        along = -along
        s = -s
    low, high = np.quantile(s, [0.01, 0.99])
    span = high - low
    if span <= 0:
        raise ValueError("Degenerate interface")
    u = (s - low) / span
    v = (raw - origin) @ across
    edges = np.linspace(0, 1, SETTINGS["bins"] + 1)
    uv = np.array(
        [
            (np.median(u[k]), np.median(v[k]))
            for a, b in zip(edges[:-1], edges[1:])
            if (k := ((u >= a) & (u < b))).any()
        ]
    )
    if len(uv) < 10:
        raise ValueError("Insufficient density-interface extent")
    knots = np.r_[np.zeros(4), [0.25, 0.5, 0.75], np.ones(4)]
    ncoef = len(knots) - 4
    basis = BSpline(knots, np.eye(ncoef), 3)
    B = basis(uv[:, 0])
    grid = np.linspace(-0.5, 1.5, 161)
    G = basis(grid)
    D2 = basis.derivative(2)(grid) / span**2
    smooth = basis.derivative(2)(np.linspace(0, 1, 81))
    L = SETTINGS["smoothing_penalty"] * (smooth.T @ smooth) / len(smooth)
    rows = features.copy()
    rows["distance_to_interface_um"] = cKDTree(raw).query(
        rows[["x_px", "y_px"]].to_numpy() * px
    )[0]
    coords = rows[["x_px", "y_px"]].to_numpy() * px - origin
    fu = (coords @ along - low) / span
    fv = coords @ across
    eligible = (
        rows.feature_valid.to_numpy()
        & (rows.brightness_weight.to_numpy() > 0)
        & (rows.distance_to_interface_um.to_numpy() <= SETTINGS["corridor_um"])
        & (fu >= 0)
        & (fu <= 1)
    )
    rows["af488_geometry_eligible"] = eligible
    rows["geometry_weight"] = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        k = eligible & (fu >= a) & (fu <= b if b == 1 else fu < b)
        w = rows.loc[k, "brightness_weight"].to_numpy()
        if len(w):
            rows.loc[k, "geometry_weight"] = (
                SETTINGS["af488_weight_per_bin"] * w / max(1.0, w.sum())
            )
    bw = rows.geometry_weight.to_numpy()
    active = bw > 0
    F = basis(fu[active])
    w = bw[active]
    target = fv[active]
    results = []
    initial = np.linalg.lstsq(B, uv[:, 1], rcond=None)[0]
    for sign in [-1, 1]:
        H = B.T @ B + L
        rhs = B.T @ uv[:, 1]

        def solve(H, rhs, extra=()):
            objective = lambda c: float(0.5 * c @ H @ c - rhs @ c)
            jac = lambda c: H @ c - rhs
            return minimize(
                objective,
                initial,
                jac=jac,
                method="SLSQP",
                constraints=[LinearConstraint(sign * D2, 0, np.inf), *extra],
                options={"maxiter": 500, "ftol": 1e-09},
            )

        base = solve(H, rhs)
        if not base.success:
            continue
        maxshift = SETTINGS["max_shift_from_density_um"]
        anchor = G @ base.x
        refined = (
            solve(
                H + F.T @ (w[:, None] * F),
                rhs + F.T @ (w * target),
                [LinearConstraint(G, anchor - maxshift, anchor + maxshift)],
            )
            if active.any()
            else base
        )
        if not refined.success:
            continue
        c = refined.x
        density_loss = float(np.mean((B @ base.x - uv[:, 1]) ** 2))
        results.append((density_loss, c, base.x, sign))
    if not results:
        raise ValueError("Constrained smooth fit did not converge")
    _, c, base, sign = min(results, key=lambda r: r[0])
    greville = np.array([knots[i + 1 : i + 4].mean() for i in range(ncoef)])
    controls = origin + (low + span * greville[:, None]) * along + c[:, None] * across
    model = BSpline(knots, controls, 3)
    residual = np.linalg.norm(model(np.clip(u, 0, 1)) - raw, axis=1)
    med, p95 = np.quantile(residual, [0.5, 0.95])
    shift = float(np.max(np.abs(G @ (c - base))))
    info = dict(
        geometry_model="convex cubic graph in density-interface PCA coordinates",
        interface_fit_median_distance_um=float(med),
        interface_fit_p95_distance_um=float(p95),
        fit_qc_passed=bool(med <= 10 and p95 <= 20),
        af488_eligible_somata=int(active.sum()),
        af488_total_weight=float(w.sum()),
        max_AF488_shift_um=shift,
        curvature_sign=int(sign),
        settings=SETTINGS,
        review_required=True,
    )
    return (model, info, rows)
