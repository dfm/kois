#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["LightCurve", "Model"]

import numpy as np
from bart.data import LightCurve
from . import _kois


class KOILightCurve(LightCurve):

    def remove_polynomial(self, periods, epochs, durations, order=1):
        t = self.time
        m = np.ones_like(t, dtype=bool)
        for p, t0, dt in zip(periods, epochs, durations):
            hp = 0.5 * p
            m[np.abs((t - t0 + hp) % p - hp) < dt] = 0

        if np.sum(m) < 10:
            return False

        p = np.polyfit(t[m], self.flux[m], order, w=self.ivar[m])
        self.flux /= np.polyval(p, t)
        return True

    def lnlike(self, light_curve):
        return np.sum(-0.5 * (light_curve - self.flux) ** 2 * self.ivar)


class Model(object):

    def __init__(self, ldp, epoch_tol=1.0, period_tol=7e-4, tol=1e-7,
                 max_depth=4):
        self.ldp = ldp
        self.datasets = []

        self.fstar = 1.0
        self.periods = np.empty(0)
        self.epochs = np.empty(0)
        self.durations = np.empty(0)
        self.rors = np.empty(0)
        self.impacts = np.empty(0)

        self.epoch_tol = epoch_tol
        self.period_tol = period_tol

        self.tol = tol
        self.max_depth = max_depth

    @property
    def vector(self):
        return np.concatenate(([self.fstar, self.ldp.q1, self.ldp.q2],
                               self.periods, self.epochs, self.durations,
                               self.rors, self.impacts))

    @vector.setter
    def vector(self, v):
        self.fstar = v[0]
        self.ldp.q1, self.ldp.q2 = v[1:3]
        self.periods, self.epochs, self.durations, self.rors, self.impacts \
            = v[3:].reshape([-1, len(self.periods)])

    def add_koi(self, period, epoch, duration, ror, impact):
        self.periods = np.append(self.periods, period)
        self.epochs = np.append(self.epochs, epoch)
        self.durations = np.append(self.durations, duration)
        self.rors = np.append(self.rors, ror)
        self.impacts = np.append(self.impacts, impact)

        # Save the initial values for use in the prior.
        self.initial_periods = np.array(self.periods)
        self.initial_epochs = np.array(self.epochs)
        self.initial_durations = np.array(self.durations)

    def get_light_curve(self, t, K=1, texp=0):
        mu1, mu2 = self.ldp.coeffs
        return _kois.light_curve(t, texp, self.periods, self.epochs,
                                 self.durations, self.rors, self.impacts,
                                 mu1, mu2, self.tol, self.max_depth)*self.fstar

    def lnlike(self):
        return sum([lc.lnlike(self.get_light_curve(lc.time, lc.K, lc.texp))
                    for lc in self.datasets])

    def lnprior(self):
        # Check the physicality of limb-darkening coefficients.
        if not (0.0 <= self.ldp.q1 <= 1.0 and 0.0 <= self.ldp.q2 <= 1.0):
            return -np.inf

        # Physical priors on other parameters.
        if not 0 < self.fstar < 2:
            return -np.inf
        if np.any(self.periods < 0):
            return -np.inf
        if (np.any(self.epochs < -self.periods) or
                np.any(self.epochs > self.periods)):
            return -np.inf
        if np.any(self.impacts < 0) or np.any(self.impacts > 2):
            return -np.inf
        if np.any(self.durations < 0) or np.any(self.durations > 10):
            return -np.inf
        if np.any(self.rors < 1e-4) or np.any(self.rors > 1):
            return -np.inf

        # Priors based on trimming of the datasets.
        if np.any(np.abs(self.epochs-self.initial_epochs) >
                  self.epoch_tol * self.initial_durations):
            return -np.inf
        if np.any(np.abs(self.periods-self.initial_periods)
                  / self.initial_periods >
                  self.period_tol * self.initial_durations):
            return -np.inf

        return -2*np.sum(np.log(self.rors))

    def lnprob(self):
        lp = self.lnprior()
        if not np.isfinite(lp):
            return -np.inf
        ll = self.lnlike()
        return ll + lp

    def __call__(self, p):
        self.vector = p
        return self.lnprob()
