import numpy as np
import scipy.signal


def estimate_velocity(paths, frame, delay):
    """
    Estimate time-average velocities by linear regression on paths.
    """

    # We'd like to estimate velocity in a symmetric window around the specified
    # time point.
    window = delay + 1
    back_window = window // 2
    forw_window = window - back_window

    if frame - back_window < 0:
        back_window = frame
    if frame + forw_window > len(paths):
        forw_window = len(paths) - frame

    window = back_window + forw_window

    # Fit x(t) = a + v t to the paths. The slope v gives the velocity.
    times = np.arange(window)
    centered_times = times - times.mean(0)
    cov_to_slope = 1 / np.sum(np.square(centered_times))

    paths = paths[frame - back_window:frame + forw_window]
    centered_paths = paths - paths.mean(0)

    velocities = np.einsum("t,tik->ik", centered_times, centered_paths) * cov_to_slope

    return velocities


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
