import numpy as np


def probability(counts: np.ndarray) -> np.ndarray:

    prob = np.zeros_like(counts)
    prob[..., 0, 0] = 1

    counts_s = counts.sum((-1, -2))[..., None, None]
    return np.divide(
        counts,
        counts_s,
        where=(counts_s > 0),
        out=prob,
    )


def conditional(counts: np.ndarray) -> np.ndarray:

    prob = probability(counts)
    joint = prob[..., 1, 1]
    marginal = prob[..., 1, 1] + prob[..., 1, 0]

    cond = np.divide(
        joint,
        marginal,
        where=(marginal > 0),
        out=np.zeros_like(joint)
    )
    return cond


def covariance(counts: np.ndarray) -> np.ndarray:

    prob = probability(counts)
    prob = prob / prob.sum((-2, -1))[..., None, None]

    ex = (prob[..., 1, 1] + prob[..., 1, 0]).mean(1)
    ey = (prob[..., 1, 1] + prob[..., 0, 1]).mean(0)
    exy = prob[..., 1, 1]

    return exy - ex[None, :] * ey[:, None]


def correlation(counts: np.ndarray) -> np.ndarray:

    cov = covariance(counts)
    std = np.sqrt(np.diag(cov))

    cov = cov / (std[:, None] * std[None, :])
    return np.nan_to_num(cov, 0.0)


def cooccurrence(counts: np.ndarray) -> np.ndarray:

    counts_s = counts.sum((-2, -1))
    return counts_s / counts_s.sum(1)
