import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from scipy.interpolate import BSpline
from scipy.optimize import minimize, LinearConstraint, Bounds


def drop_width(u, p, side):
    v = u * side
    ok = np.isfinite(p) & (v >= 0) & (v <= 45)
    v = v[ok]
    q = p[ok]
    o = np.argsort(v)
    v = v[o]
    q = q[o]
    if len(v) < 15 or v[0] > 1:
        return np.nan
    q = ndi.gaussian_filter1d(q, 1.5)
    k = np.where(v <= 8)[0][np.argmax(q[v <= 8])]
    peak = q[k]
    if peak <= 1e-12:
        return np.nan
    for j in range(k + 1, len(v) - 4):
        if 5 <= v[j] <= 40 and np.all(q[j : j + 5] <= 0.6 * peak):
            return v[j]
    return np.nan


def smooth_width(s, raw, lower=None):
    valid = np.isfinite(raw)
    target = (
        np.interp(s, s[valid], raw[valid])
        if valid.sum() >= 3
        else np.full(len(s), 15.0)
    )
    w = np.where(valid, 1.0, 0.1)
    knots = np.r_[
        np.repeat(s[0], 4), np.arange(s[0] + 60, s[-1], 60), np.repeat(s[-1], 4)
    ]
    n = len(knots) - 4
    basis = BSpline(knots, np.eye(n), 3)
    B = basis(s)
    grid = np.linspace(s[0], s[-1], max(2, len(s) * 8))
    G = basis(grid)
    D = basis.derivative(1)(grid)
    D2 = basis.derivative(2)(grid)
    constraints = [LinearConstraint(D, -0.35, 0.35)]
    if lower is not None:
        constraints.append(LinearConstraint(B, lower, np.inf))
        start = np.full(n, min(60.0, max(5.0, np.max(lower) + 0.1)))
    else:
        start = np.full(n, 15.0)

    def fun(c):
        e = (B @ c - target) / 3.0
        k = D2 @ c
        return np.mean(w * 9 * (np.sqrt(1 + e * e) - 1)) + 25.0**4 * np.mean(k * k)

    def jac(c):
        e = (B @ c - target) / 3.0
        return B.T @ (w * 3 * e / np.sqrt(1 + e * e)) / len(
            s
        ) + 2 * 25.0**4 * D2.T @ (D2 @ c) / len(grid)

    fit = minimize(
        fun,
        start,
        jac=jac,
        method="SLSQP",
        bounds=Bounds(5.0, 60.0),
        constraints=constraints,
        options={"maxiter": 1000, "ftol": 1e-09},
    )
    if not fit.success:
        raise RuntimeError("Smooth width solver failed: " + fit.message)
    values = B @ fit.x
    if lower is not None and np.max(lower - values) > 1e-05:
        raise RuntimeError("Soma envelope constraint failed")
    return (
        values,
        {
            "valid_profiles": int(valid.sum()),
            "default_used": bool(valid.sum() < 3),
            "coefficients": fit.x.tolist(),
            "knots_um": knots.tolist(),
        },
    )


def signed_projection(points, c, normal, tree):
    _, i = tree.query(points)
    return (i, np.sum((points - c[i]) * normal[i], axis=1))


def compute_band(soma, nuc, af, dapi, px, white, rows):
    curve = white * px
    arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(curve, axis=0), axis=1))]
    s = np.arange(0, arc[-1], 4.0)
    c = np.column_stack([np.interp(s, arc, curve[:, j]) for j in [0, 1]])
    t = np.gradient(c, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    normal = np.column_stack([-t[:, 1], t[:, 0]])
    ds = 4
    shape = ((nuc.shape[0] + 3) // 4, (nuc.shape[1] + 3) // 4)
    density = np.zeros(shape)
    ids = np.unique(nuc)
    ids = ids[ids != 0]
    centers = ndi.center_of_mass(np.ones(nuc.shape, np.uint8), nuc, ids)
    for y, x in centers:
        density[min(int(y / ds), shape[0] - 1), min(int(x / ds), shape[1] - 1)] += 1
    density = ndi.gaussian_filter(density, 6 / (px * ds))
    u = np.arange(-80, 81.0)

    def sample():
        xy = c[:, None, :] + normal[:, None, :] * u[None, :, None]
        return ndi.map_coordinates(
            density,
            [xy[:, :, 1] / (px * ds), xy[:, :, 0] / (px * ds)],
            order=1,
            mode="constant",
            cval=np.nan,
        )

    prof = sample()
    if np.nanmedian(prof[:, u == 30] - prof[:, u == -30]) < 0:
        normal = -normal
        prof = sample()
    # Median of neighboring parallel profiles suppresses single-nucleus gaps.
    consensus = ndi.generic_filter(
        prof,
        lambda v: np.median(v[np.isfinite(v)]) if np.isfinite(v).any() else np.nan,
        size=(7, 1),
        mode="nearest",
    )
    raw_s = np.array([drop_width(u, p, -1) for p in consensus])
    raw_d = np.array([drop_width(u, p, 1) for p in consensus])
    core_s, sm = smooth_width(s, raw_s)
    core_d, dm = smooth_width(s, raw_d)
    border_ids = set(
        np.unique(np.concatenate([soma[0], soma[-1], soma[:, 0], soma[:, -1]]))
    ) - {0}
    tree = cKDTree(c)
    boxes = ndi.find_objects(soma)
    lower = core_s.copy()
    support = []
    review = []
    # Exploratory soma-support cohort, not a cell-type classifier.
    for r in rows:
        area = float(r["soma_area_um2"])
        z = float(r["cytoplasmic_contrast"])
        sid = int(r["soma_id"])
        if int(r["soma_id"]) in border_ids:
            continue
        if not (area >= 65 and np.isfinite(z) and z >= 3):
            continue
        box = boxes[sid - 1]
        yy, xx = np.where(soma[box] == sid)
        points = np.column_stack([xx + box[1].start, yy + box[0].start]) * px
        i, d = signed_projection(points, c, normal, tree)
        core = (d >= -core_s[i]) & (d <= core_d[i])
        fraction = float(np.mean(core))
        gap = np.maximum.reduce([-core_s[i] - d, d - core_d[i], np.zeros(len(d))])
        nearest_gap = float(np.min(gap))
        eligible = fraction >= 0.1 or nearest_gap <= 5.0
        record = dict(
            soma_id=sid,
            area_um2=area,
            cytoplasmic_contrast=z,
            overlap_density_core=fraction,
            nucleus_assignment=r["primary_nucleus_status"],
            support_eligible=eligible,
            nearest_gap_to_density_core_um=nearest_gap,
        )
        if eligible:
            for k in np.unique(i):
                lower[k] = max(lower[k], float(np.max(-d[i == k])) + 1.0)
            support.append(record)
        else:
            review.append(record)
    if lower.max() > 59.9:
        raise RuntimeError(
            "Required extension exceeds 60 um; review support cells rather than silently widening"
        )
    roi_s, rm = smooth_width(s, core_s, lower)
    # Entire masks, not only the portion inside a band, are the unit for later puncta analysis.
    records = []
    coremask = np.zeros(soma.shape, bool)
    roimask = np.zeros(soma.shape, bool)
    xx = np.arange(soma.shape[1]) * px
    for y0 in range(0, soma.shape[0], 128):
        y1 = min(y0 + 128, soma.shape[0])
        yy = np.arange(y0, y1) * px
        X, Y = np.meshgrid(xx, yy)
        pts = np.column_stack([X.ravel(), Y.ravel()])
        i, d = signed_projection(pts, c, normal, tree)
        coremask[y0:y1] = ((d >= -core_s[i]) & (d <= core_d[i])).reshape(y1 - y0, -1)
        roimask[y0:y1] = ((d >= -roi_s[i]) & (d <= core_d[i])).reshape(y1 - y0, -1)
    supportids = {r["soma_id"] for r in support}
    for r in rows:
        sid = int(r["soma_id"])
        b = boxes[sid - 1]
        m = soma[b] == sid
        r = dict(r)
        r.update(
            density_core_fraction=float(np.mean(coremask[b][m])),
            analysis_ROI_fraction=float(np.mean(roimask[b][m])),
            soma_support=sid in supportids,
        )
        records.append(r)
    return (
        roimask,
        coremask,
        records,
        np.column_stack([s, c / px, normal, core_s, core_d, roi_s, core_d]),
        {"sparse_width_fit": sm, "dense_width_fit": dm, "soma_envelope_fit": rm},
    )
