#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["load_system", "LightCurve", "Model"]

import numpy as np

import kplr
from kplr.ld import get_quad_coeffs

import bart
from bart.data import LightCurve

from . import _kois


def load_system(koi_id):
    client = kplr.API()
    kois = client.koi("{0}.01".format(koi_id)).star.kois

    # Set up the initial limb-darkening profile.
    mu1, mu2 = get_quad_coeffs(kois[0].koi_steff)
    ldp = bart.ld.QuadraticLimbDarkening(mu1, mu2, 100)

    # Set up the model.
    model = Model(ldp)

    # Loop over KOIs and add the initial values to the model.
    for koi in kois:
        P = koi.koi_period

        # Skip single transits.
        if P < 0:
            continue

        t0 = koi.koi_time0bk % P
        model.add_koi(P, t0, koi.koi_duration/24., koi.koi_ror, koi.koi_impact)

    # Load the light curves.
    lcs = koi.get_light_curves(short_cadence=False)
    datasets = []
    for lc in lcs:
        data = lc.read()
        m = data["SAP_QUALITY"] == 0
        data = KOILightCurve(data["TIME"][m],
                             data["SAP_FLUX"][m],
                             data["SAP_FLUX_ERR"][m])
        datasets += data.active_window(model.periods, model.epochs,
                                       4*model.durations).autosplit()

    # Remove datasets with not enough data points.
    datasets = [d for d in datasets if len(d.time) > 50]

    # De-trend the data.
    model.datasets += [d for d in datasets
                       if d.remove_polynomial(model.periods, model.epochs,
                                              2*model.durations, 1)]

    return model


class KOILightCurve(LightCurve):

    def remove_polynomial(self, periods, epochs, durations, order=1):
        t = self.time
        m = np.ones_like(t, dtype=bool)
        for p, t0, dt in zip(periods, epochs, durations):
            hp = 0.5 * p
            m[np.abs((t - t0 + hp) % p - hp) < dt] = 0

        if not np.sum(m):
            return False

        p = np.polyfit(t[m], self.flux[m], order, w=self.ivar[m])
        self.flux /= np.polyval(p, t)
        return True

    def lnlike(self, light_curve):
        return np.sum(-0.5 * (light_curve - self.flux) ** 2 * self.ivar)


class Model(object):

    def __init__(self, ldp):
        self.ldp = ldp
        self.datasets = []

        self.periods = np.empty(0)
        self.epochs = np.empty(0)
        self.durations = np.empty(0)
        self.rors = np.empty(0)
        self.impacts = np.empty(0)

    @property
    def vector(self):
        return np.concatenate(([self.ldp.gamma1, self.ldp.gamma2],
                               self.periods, self.epochs, self.durations,
                               self.rors, self.impacts))

    @vector.setter
    def vector(self, v):
        p = np.array(self.periods)
        self.ldp.gamma1, self.ldp.gamma2 = v[:2]
        self.periods, self.epochs, self.durations, self.rors, self.impacts \
            = v[2:].reshape([-1, len(self.periods)])

    def add_koi(self, period, epoch, duration, ror, impact):
        self.periods = np.append(self.periods, period)
        self.epochs = np.append(self.epochs, epoch)
        self.durations = np.append(self.durations, duration)
        self.rors = np.append(self.rors, ror)
        self.impacts = np.append(self.impacts, impact)

    def get_light_curve(self, t, K=1, texp=0):
        return _kois.light_curve(t, K, texp, self.periods, self.epochs,
                                 self.durations, self.rors, self.impacts,
                                 self.ldp.bins, self.ldp.intensity)

    def lnlike(self):
        return sum([lc.lnlike(self.get_light_curve(lc.time, lc.K, lc.texp))
                    for lc in self.datasets])

    def lnprior(self):
        u1, u2 = self.ldp.gamma1, self.ldp.gamma2
        if u1 < 0.0 or u1+u2 > 1.0 or u1+2*u2 < 0:
            return -np.inf

        if np.any(self.impacts < 0) or np.any(self.impacts > 1):
            return -np.inf

        if np.any(self.durations < 0):
            return -np.inf

        return 0.0

    def lnprob(self):
        lp = self.lnprior()
        if not np.isfinite(lp):
            return -np.inf
        ll = self.lnlike()
        return ll + lp

    def __call__(self, p):
        self.vector = p
        return self.lnprob()
