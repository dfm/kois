#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["load_system", "LightCurve", "Model"]

import logging
import numpy as np

import kplr
from kplr.ld import get_quad_coeffs

import bart
from bart.data import LightCurve

from . import _kois


def load_system(kepid, lc_window_factor=4, sc_window_factor=4,
                detrend_window_factor=1.5, min_dataset_size=10, poly_order=1,
                ldp_nbins=100):
    client = kplr.API()

    logging.info("Getting the catalog parameters for KIC {0}"
                 .format(kepid))
    kois = sorted([k for k in client.star(kepid).kois
                   if k.koi_disposition in ["CANDIDATE", "CONFIRMED"]],
                  key=lambda k: k.kepoi_name)
    if not len(kois):
        raise RuntimeError("No KOIs listed in the catalog for ID '{0}'"
                           .format(kepid))

    # Set up the initial limb-darkening profile.
    teff = kois[0].koi_steff
    mu1, mu2 = get_quad_coeffs(teff if teff is not None else 5778)
    mu1 = max(0.0, mu1)
    ldp = bart.ld.QuadraticLimbDarkening(mu1, mu2, ldp_nbins)

    # Set up the model.
    model = Model(ldp, epoch_tol=sc_window_factor,
                  period_tol=sc_window_factor*3e-4)

    # Loop over KOIs and add the initial values to the model.
    for koi in kois:
        # Parse the KOI parameters.
        P = koi.koi_period
        t0 = koi.koi_time0bk % P
        ror = koi.koi_ror
        b = koi.koi_impact

        # Skip single transits.
        if P < 0:
            continue

        # Convert the KOI duration to the b=0 duration that we need.
        tau = (koi.koi_duration/24.) / np.sqrt((1+ror)**2 - b*b)

        # Add the KOI to the model.
        model.add_koi(P, t0, tau, ror, b)

    # Load the light curves.
    logging.info("Loading datasets")
    lcs = koi.get_light_curves()
    datasets = []
    for lc in lcs:
        data = lc.read()
        m = data["SAP_QUALITY"] == 0
        texp = (kplr.EXPOSURE_TIMES[1] if lc.sci_archive_class == "CLC"
                else kplr.EXPOSURE_TIMES[0])
        factor = (lc_window_factor if lc.sci_archive_class == "CLC"
                  else sc_window_factor)
        data = KOILightCurve(data["TIME"][m],
                             data["PDCSAP_FLUX"][m],
                             data["PDCSAP_FLUX_ERR"][m], texp=texp,
                             K=5 if lc.sci_archive_class == "CLC" else 3)
        datasets += (data.active_window(model.periods, model.epochs,
                                        factor*model.durations)
                     .autosplit())

    # Remove datasets with not enough data points.
    datasets = [d for d in datasets if len(d.time) > min_dataset_size]
    if not len(datasets):
        raise RuntimeError("No datasets survived the cuts.")
    logging.info("{0} datasets were found and trimmed".format(len(datasets)))

    # De-trend the data.
    model.datasets += [d for d in datasets
                       if d.remove_polynomial(model.periods, model.epochs,
                                              detrend_window_factor
                                              * model.durations, poly_order)]
    logging.info("{0} datasets were de-trended".format(len(model.datasets)))
    if not len(datasets):
        raise RuntimeError("No datasets could be de-trended.")

    return kois[0].kepoi_name.split(".")[0], model


class KOILightCurve(LightCurve):

    def remove_polynomial(self, periods, epochs, durations, order=1):
        t = self.time
        m = np.ones_like(t, dtype=bool)
        for p, t0, dt in zip(periods, epochs, durations):
            hp = 0.5 * p
            m[np.abs((t - t0 + hp) % p - hp) < dt] = 0

        if not np.sum(m):
            return False

        p = np.polyfit(t[m], self.flux[m], order)  # , w=self.ivar[m])
        self.flux /= np.polyval(p, t)
        return True

    def lnlike(self, light_curve):
        return np.sum(-0.5 * (light_curve - self.flux) ** 2 * self.ivar)


class Model(object):

    def __init__(self, ldp, epoch_tol=1.0, period_tol=7e-4):
        self.ldp = ldp
        self.datasets = []

        self.periods = np.empty(0)
        self.epochs = np.empty(0)
        self.durations = np.empty(0)
        self.rors = np.empty(0)
        self.impacts = np.empty(0)

        self.epoch_tol = epoch_tol
        self.period_tol = period_tol

    @property
    def vector(self):
        return np.concatenate(([self.ldp.gamma1, self.ldp.gamma2],
                               self.periods, self.epochs, self.durations,
                               self.rors, self.impacts))

    @vector.setter
    def vector(self, v):
        self.ldp.gamma1, self.ldp.gamma2 = v[:2]
        self.periods, self.epochs, self.durations, self.rors, self.impacts \
            = v[2:].reshape([-1, len(self.periods)])

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
        return _kois.light_curve(t, K, texp, self.periods, self.epochs,
                                 self.durations, self.rors, self.impacts,
                                 self.ldp.bins, self.ldp.intensity)

    def lnlike(self):
        return sum([lc.lnlike(self.get_light_curve(lc.time, lc.K, lc.texp))
                    for lc in self.datasets])

    def lnprior(self):
        # Physical limb-darkening priors.
        u1, u2 = self.ldp.gamma1, self.ldp.gamma2
        if u1 < 0.0 or u1+u2 > 1.0 or u1+2*u2 < 0:
            return -np.inf

        # Physical priors on other parameters.
        if np.any(self.periods < 0):
            return -np.inf
        if (np.any(self.epochs < -self.periods) or
                np.any(self.epochs > self.periods)):
            return -np.inf
        if np.any(self.impacts < 0) or np.any(self.impacts-self.rors > 1):
            return -np.inf
        if np.any(self.durations < 0):
            return -np.inf
        if np.any(self.rors < 0) or np.any(self.rors > 1):
            return -np.inf

        # Priors based on trimming of the datasets.
        if np.any(np.abs(self.epochs-self.initial_epochs) >
                  self.epoch_tol * self.initial_durations):
            return -np.inf
        if np.any(np.abs(self.periods-self.initial_periods)
                  / self.initial_periods >
                  self.period_tol * self.initial_durations):
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
