import numpy as np
from sklearn.mixture import GaussianMixture


def fit_model(X, features, seed=17):
    info = {"n": len(X), "supported": False, "features": features}
    if len(X) < 30:
        return (None, info)
    med = np.median(X, axis=0)
    scale = 1.4826 * np.median(abs(X - med), axis=0)
    scale = np.where(scale > 1e-06, scale, np.std(X, axis=0))
    scale = np.where(scale > 1e-06, scale, 1)
    Z = (X - med) / scale
    a = GaussianMixture(1, reg_covar=0.02, n_init=10, random_state=seed).fit(Z)
    b = GaussianMixture(2, reg_covar=0.02, n_init=10, random_state=seed).fit(Z)
    size = [
        features.index(k)
        for k in ["log_soma_area", "log_nucleus_area"]
        if k in features
    ]
    target = int(np.argmax(np.mean(b.means_[:, size], axis=1)))
    other = 1 - target
    counts = np.bincount(b.predict(Z), minlength=2)
    delta = float(a.bic(Z) - b.bic(Z))
    direction = bool(np.all(b.means_[target, size] > b.means_[other, size]))
    if "log1p_AF488" in features:
        direction = direction and bool(
            b.means_[target, features.index("log1p_AF488")]
            > b.means_[other, features.index("log1p_AF488")]
        )
    threshold = float(
        np.exp((b.means_[target, 0] + b.means_[other, 0]) / 2 * scale[0] + med[0])
    )
    info.update(
        supported=bool(delta >= 10 and counts.min() >= 10 and direction),
        BIC_gain=delta,
        component_counts=counts.tolist(),
        soma_area_midpoint_um2=threshold,
        target_component=target,
        means_standardized=b.means_.tolist(),
        median=med.tolist(),
        scale=scale.tolist(),
    )
    return ((b, med, scale, target, threshold), info)


def predict(X, model):
    b, med, scale, target, threshold = model
    p = b.predict_proba((X - med) / scale)[:, target]
    positive = (p >= 0.8) & (np.exp(X[:, 0]) >= threshold)
    return (p, positive)
