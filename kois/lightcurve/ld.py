#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["QuadraticLimbDarkening"]

import numpy as np


class QuadraticLimbDarkening(object):
    """
    Based on code from bart.

    """

    def __init__(self, gamma1, gamma2, bins=500):
        self.coeffs = (gamma1, gamma2)
        try:
            len(bins)
        except TypeError:
            bins = np.sqrt(np.linspace(0, 1, int(bins)+1))[1:]
        assert np.all(bins > 0) and np.all(bins <= 1)
        self.bins = np.sort(bins)
        self.bins /= self.bins[-1]
        self._intensity = None

    @property
    def q1(self):
        return self._q1

    @q1.setter
    def q1(self, v):
        self._intensity = None
        self._q1 = v

    @property
    def q2(self):
        return self._q2

    @q2.setter
    def q2(self, v):
        self._intensity = None
        self._q2 = v

    @property
    def coeffs(self):
        q1, q2 = self._q1, self._q2
        q1 = np.sqrt(np.abs(q1))
        return 2*q1*q2, q1*(1-2*q2)

    @coeffs.setter
    def coeffs(self, value):
        u1, u2 = value
        u2 = u1+u2
        self.q1, self.q2 = u2*u2, 0.5*u1/u2

    @property
    def intensity(self):
        if self._intensity is None:
            b = np.append(0, self.bins)
            self._intensity = self.integrate(b[:-1], b[1:])
            self._intensity /= np.pi * (b[1:]**2 - b[:-1]**2)
        return self._intensity

    def integrate(self, a, b):
        m1, m2 = self.coeffs
        a2, b2 = a * a, b * b
        th = 0.5 * 3
        k1 = 0.5 * (b2 - a2) * (1 - m1 - 2 * m2)
        k2 = (m1 + 2 * m2) * ((1 - a2) ** th - (1 - b2) ** th) / 3
        k3 = 0.25 * m2 * (b2 * b2 - a2 * a2)
        return 2 * np.pi * (k1 + k2 + k3)

    def __call__(self, r):
        onemmu = 1 - np.sqrt(1 - r * r)
        m1, m2 = self.coeffs
        return 1 - m1 * onemmu - m2 * onemmu * onemmu


if __name__ == "__main__":
    import matplotlib.pyplot as pl

    c = np.array([0.3, 0.1])
    ld = QuadraticLimbDarkening(*c)
    print(ld.coeffs - c)
    pl.plot(ld.bins, ld.intensity)
    pl.plot(ld.bins, ld(ld.bins))
    pl.savefig("test.png")
