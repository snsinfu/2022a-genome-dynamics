import numpy as np
import scipy.interpolate


def get_refinement_method(name):
    return refinement_methods[name]


def refine_path_spline(path, n):
    """
    Refine path to n points using cubic spline interpolation.
    """
    # This parameterization ensures that each sample point is the midpoint of
    # a bin.
    u = (np.arange(len(path)) + 0.5) / len(path)
    fine_u = (np.arange(n) + 0.5) / n
    tck, _ = scipy.interpolate.splprep(path.T, u=u, k=3, s=0)
    fine_path = np.transpose(scipy.interpolate.splev(fine_u, tck))
    return fine_path


refinement_methods = {
    "spline": refine_path_spline,
}
