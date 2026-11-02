import numpy as np
from scipy.stats import multivariate_normal
from scipy.special import logsumexp


def marginal_membership(area, contrast, model):
    b, med, scale, target, _ = model
    cols = [0, 2]
    v = (np.array([np.log(area), np.log1p(contrast)]) - med[cols]) / scale[cols]
    ll = np.array(
        [
            np.log(b.weights_[k])
            + multivariate_normal.logpdf(
                v, mean=b.means_[k, cols], cov=b.covariances_[k][np.ix_(cols, cols)]
            )
            for k in range(2)
        ]
    )
    return float(np.exp(ll[target] - logsumexp(ll)))
