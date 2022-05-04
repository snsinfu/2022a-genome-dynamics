import numpy as np
import scipy.signal


def gaussian_smooth(signal, window, sigma=3):
    """
    Apply Gaussian FIR filter to signal.
    """
    rank = len(signal.shape)

    kernel = np.exp(-np.linspace(-sigma, sigma, num=window)**2 / 2)
    kernel = kernel / np.sum(kernel)
    kernel = kernel.reshape(kernel.shape + (1, ) * (rank - 1))

    lpad = window // 2
    rpad = window - lpad - 1
    pad_sizes = [(lpad, rpad)] + [(0, 0)] * (rank - 1)
    padded_signal = np.pad(signal, pad_sizes, mode="reflect")

    return scipy.signal.fftconvolve(padded_signal, kernel, mode="valid", axes=0)
