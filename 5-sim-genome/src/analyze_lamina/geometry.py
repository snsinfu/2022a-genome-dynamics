import numpy as np


EPSILON = 1e-6


class Ellipsoid:
    def __init__(self, semiaxes):
        semiaxes = np.array(semiaxes, dtype=np.float)
        self._semiaxes = semiaxes
        self._semiaxes_invsq = semiaxes ** -2

    def distance_from_surface(self, points):
        """
        Given array of points, return second-order approximate distances from
        the surface.
        """
        s1 = self._semiaxes_invsq[None, :] * points
        s2 = self._semiaxes_invsq[None, :] * s1
        s3 = self._semiaxes_invsq[None, :] * s2

        a = np.sum(s3 * points, axis=1)
        b = np.sum(s2 * points, axis=1)
        c = np.sum(s1 * points, axis=1) - 1
        u = (b - np.sqrt(b ** 2 - a * c)) / (a + EPSILON)
        v = np.linalg.norm(s1, axis=1)

        return np.abs(u * v)
