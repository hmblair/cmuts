import numpy as np


def _validate_shape(data: np.ndarray) -> np.ndarray:

    return (
        data.shape[0] == data.shape[1] and
        data.shape[2] == 2 and
        data.shape[3] == 2
    )


def zero_diag(vals: np.ndarray, n: int = 1):

    ix = range(vals.shape[0])

    if n > 0:
        vals[ix, ix] = 0.0

    for jx in range(1, n):
        vals[ix[jx:], ix[:-jx]] = 0.0
        vals[ix[:-jx], ix[jx:]] = 0.0

    return vals


def probability(counts: np.ndarray) -> np.ndarray:

    if not _validate_shape(counts):
        raise ValueError("The data has an incorrect shape.")

    prob = np.zeros_like(counts)
    prob[..., 0, 0] = 1

    counts_s = counts.sum((-1, -2))[..., None, None]

    return np.divide(
        counts,
        counts_s,
        where=(counts_s > 0),
        out=prob,
    )


def conditional(mod: np.ndarray) -> np.ndarray:

    prob = probability(mod)

    joint = prob[..., 1, 1]
    marginal = prob[..., 1, 1] + prob[..., 1, 0]

    return np.divide(
        joint,
        marginal,
        where=(marginal > 0),
        out=np.zeros_like(joint)
    )


def covariance(mod: np.ndarray) -> np.ndarray:

    prob = probability(mod)
    prob = prob / prob.sum((-2, -1))[..., None, None]

    ex = (prob[..., 1, 1] + prob[..., 1, 0]).mean(1)
    ey = (prob[..., 1, 1] + prob[..., 0, 1]).mean(0)
    exy = prob[..., 1, 1]

    return exy - ex[None, :] * ey[:, None]


def correlation(mod: np.ndarray, nomod: np.ndarray | None = None) -> np.ndarray:

    cov = covariance(mod)
    std = np.sqrt(np.diag(cov))

    if nomod is not None:
        nmcov = covariance(nomod)
        nmstd = np.sqrt(np.diag(nmcov))
        cov = cov - nmcov
        std = np.sqrt(std ** 2 + nmstd ** 2)

    cov = cov / (std[:, None] * std[None, :])
    return np.nan_to_num(cov, 0.0)


def mutual_information(mod: np.ndarray) -> np.ndarray:

    prob = probability(mod)

    p_00 = prob[..., 0, 0]
    p_01 = prob[..., 0, 1]
    p_10 = prob[..., 1, 0]
    p_11 = prob[..., 1, 1]

    p_x0 = p_00 + p_01
    p_x1 = p_10 + p_11
    p_y0 = p_00 + p_10
    p_y1 = p_01 + p_11

    mi = np.zeros_like(p_00)

    # p(0,0) * log(p(0,0) / (p(0) * p(0)))
    valid = (p_00 > 0) & (p_x0 > 0) & (p_y0 > 0)
    mi[valid] += p_00[valid] * np.log(p_00[valid] / (p_x0[valid] * p_y0[valid]))

    # p(0,1) * log(p(0,1) / (p(0) * p(1)))
    valid = (p_01 > 0) & (p_x0 > 0) & (p_y1 > 0)
    mi[valid] += p_01[valid] * np.log(p_01[valid] / (p_x0[valid] * p_y1[valid]))

    # p(1,0) * log(p(1,0) / (p(1) * p(0)))
    valid = (p_10 > 0) & (p_x1 > 0) & (p_y0 > 0)
    mi[valid] += p_10[valid] * np.log(p_10[valid] / (p_x1[valid] * p_y0[valid]))

    # p(1,1) * log(p(1,1) / (p(1) * p(1)))
    valid = (p_11 > 0) & (p_x1 > 0) & (p_y1 > 0)
    mi[valid] += p_11[valid] * np.log(p_11[valid] / (p_x1[valid] * p_y1[valid]))

    return mi


def cooccurrence(counts: np.ndarray) -> np.ndarray:

    values = counts.sum((-2, -1))
    max = values[~np.isnan(values)].max()
    return values / max
