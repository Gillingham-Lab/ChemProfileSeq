import numpy as np
from typing import Union

def normalize(data: np.ndarray) -> np.ndarray:
    """ Normalizes the stop/del/ins/mutation data and returns a copy."""
    data = data.copy()

    # Make sure that we never divide by 0
    data[:, 0] += 1

    # Stops are normalized by shifting by 1 (0..9 / 1..10)
    data[:-1, 1] /= data[1:,0]
    # Set the not-normalized base to 0
    data[-1, 1] = 0

    # Normalize the others normally
    data[:, 2:5] /= data[:, 0, np.newaxis]

    # Set the beginning of the gene (or the end of the cDNA!) to 0
    data[:16,1:5] = 0

    return data

def weigh_for_average_reads(
        data: np.ndarray,
        treated: np.ndarray) -> np.ndarray:
    """ Weights a matrix according to the original data """
    weighted = data * treated.sum(axis=1)[:,np.newaxis]
    return weighted


def subtract(
        sample: np.ndarray,
        control: np.ndarray) -> np.ndarray:
    """ Subtracts control from sample and sets everything below 0 to 0. """
    corrected = sample - control
    corrected[corrected < 0] = 0

    return corrected


def percentile_threshold(
        input: np.ndarray,
        percentile: int = 80,
        set_to: Union[int,float] = 0,
        axis: int = 0
):
    input = input*1

    input[input < np.percentile(input, percentile, axis=axis)] = set_to
    return input


def normalize_special(theta, upper=95, lower=90):
    upper_threshold = theta < np.percentile(theta, upper)
    lower_treshold = theta > np.percentile(theta, lower)

    treshold = np.all([upper_threshold, lower_treshold], axis=0)

    theta_norm = theta / theta[treshold].mean()
    theta_norm[theta_norm > 1] = 1

    return theta_norm

def remove_common(
        sample: np.ndarray,
        control: np.ndarray
):
    result = sample * 1
    result[control > 0] = 0
    return result

def calculate_reactivities(sample, control, control_bg):
    sample_ = normalize_meanVigile2(sample)
    control_ = normalize_meanVigile2(control)

    f = 32
    R = (sample_ - 0.25 * control_) / control_bg
    R[:f] = 0
    R[-f:] = 0

    above = np.percentile(R, 95)
    below = np.percentile(R, 5)
    maximum = R[R < above].max()
    minimum = R[R > below].min()
    R -= minimum
    R /= maximum
    R[R > 1] = 1
    R[R < 0] = 0

    return R

def normalize_meanVigile2(input):
    f = 32
    mean = input[input > np.percentile(input[f:-f], 95, axis=0)].mean()
    return (input/mean)



def calculate_beta(sample, control):
    #beta = np.zeros(len(sample), dtype=np.float32)

    y = control[::-1]
    x = sample[::-1]

    y = y / control.cumsum()[::-1]
    x = x / sample.cumsum()[::-1]

    #y = control / control.sum()
    #x = sample / sample.sum()

    beta = (x - y) / (1 - y)
    #beta[beta > np.percentile(beta, 90)] = 1
    beta[beta < 0] = 0
    return beta


def weigh(
        sample: np.ndarray,
        use: int = 0,
        weighs=1.5,
        neg_weighs=1.0
):
    if isinstance(weighs, float) == False:
        assert len(weighs) == (sample.shape[1] - 1)
    else:
        weighs = [weighs] * (sample.shape[1] - 1)

    if isinstance(neg_weighs, float) == False:
        assert len(neg_weighs) == (sample.shape[1] - 1)
    else:
        neg_weighs = [neg_weighs] * (sample.shape[1] - 1)

    use_data = sample[:, use] * 1
    weigh_indecies = [x for x in range(sample.shape[1]) if x != use]

    j = 0
    for i in weigh_indecies:
        mask = sample[:, i] * 1
        mask[mask > 0] = weighs[j]
        mask[mask == 0] = neg_weighs[j]
        j += 1
        use_data *= mask

    return use_data

def _normalize(data: np.ndarray) -> np.ndarray:
    sum = data.sum(axis=0)
    sum[sum == 0] = 1
    data = data / sum
    return data