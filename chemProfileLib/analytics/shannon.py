import numpy as np

def calc_entropy(data: np.ndarray) -> np.ndarray:
    data = data + 1
    entropy = data / data.sum(axis=1)[:,np.newaxis]
    R = 2 + (entropy * np.log2(entropy)).sum(axis=1)
    entropy *= R[:, np.newaxis]

    return entropy