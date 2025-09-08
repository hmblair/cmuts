import numpy as np


def covariance(joint: np.ndarray) -> np.ndarray:

    joint = joint / joint.sum((-2, -1))[..., None, None]

    ex = (joint[..., 1, 1] + joint[..., 1, 0]).mean(1)
    ey = (joint[..., 1, 1] + joint[..., 0, 1]).mean(0)
    exy = joint[..., 1, 1]

    return exy - ex[None, :] * ey[:, None]


def correlation(joint: np.ndarray) -> np.ndarray:

    cov = covariance(joint)
    std = np.sqrt(np.diag(cov))

    cov = cov / (std[:, None] * std[None, :])
    return np.nan_to_num(cov, 0.0)
