import numpy as np
from scipy.interpolate import LSQUnivariateSpline, BSpline


def smooth_arc(raw_um, interior_knots=1, n_output=1600):
    u, r = arc_resample(raw_um)
    knots = np.r_[np.zeros(4), np.linspace(0, 1, interior_knots + 2)[1:-1], np.ones(4)]
    model = make_lsq_spline(u, r, knots, k=3)
    t = np.linspace(0, 1, n_output)
    curve = model(t)
    d1 = model.derivative(1)(t)
    d2 = model.derivative(2)(t)
    k = (d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]) / np.linalg.norm(d1, axis=1) ** 3
    return (curve, k, model)


def curvature_sign_changes(k):
    signs = np.sign(k[np.abs(k) > 1e-05])
    return int(np.count_nonzero(np.diff(signs))) if len(signs) > 1 else 0
