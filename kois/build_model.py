#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["build_model"]

import logging
import numpy as np

import kplr
from kplr.ld import get_quad_coeffs

from .lightcurve.ld import QuadraticLimbDarkening
from .kois import Model, KOILightCurve


def build_model(koi_id, window_factor=4.0, detrend_factor=1.5, poly_order=1):
    # Set up the model.
    client = kplr.API()
    koi = client.koi(koi_id)
    assert koi.koi_disposition in ["CANDIDATE", "CONFIRMED"], \
        "That KOI is not a candidate or confirmed system"

    # Limb darkening and initial model.
    teff = koi.koi_steff
    mu1, mu2 = get_quad_coeffs(teff if teff is not None else 5778)
    ldp = QuadraticLimbDarkening(mu1, mu2)
    ldp.q1 = min(0.9, max(0.1, ldp.q1))
    ldp.q2 = min(0.9, max(0.1, ldp.q2))
    model = Model("KOI {0}".format(koi_id), ldp, epoch_tol=window_factor,
                  period_tol=window_factor*3e-4)

    # Download and trim the datasets.
    model.add_koi(koi.koi_period, koi.koi_time0bk % koi.koi_period,
                  koi.koi_duration / 24.0, koi.koi_ror, koi.koi_impact)

    # Load the light curves.
    logging.info("Loading datasets")
    datasets = []
    for lc in koi.get_light_curves():
        try:
            data = lc.read()
        except:
            continue
        m = data["SAP_QUALITY"] == 0
        is_lc = lc.sci_archive_class == "CLC"
        texp = kplr.EXPOSURE_TIMES[1] if is_lc else kplr.EXPOSURE_TIMES[0]
        window = np.array([max(window_factor*d,
                           window_factor*d+5*kplr.EXPOSURE_TIMES[1]/86400.)
                           for d in model.durations])
        data = KOILightCurve(data["TIME"][m],
                             data["SAP_FLUX"][m],
                             data["SAP_FLUX_ERR"][m], texp=texp)
        data.is_lc = is_lc
        datasets += (data.active_window(model.periods, model.epochs,
                                        window)
                     .autosplit(ttol=1.5))

    # Remove datasets with not enough data points.
    datasets = [d for d in datasets if len(d.time) > 10]
    if not len(datasets):
        raise RuntimeError("No datasets survived the cuts.")
    logging.info("{0} datasets were found and trimmed".format(len(datasets)))

    # De-trend the data.
    model.datasets += [d for d in datasets
                       if d.remove_polynomial(model.periods, model.epochs,
                                              detrend_factor
                                              * model.durations,
                                              poly_order,
                                              nmin=10 if d.is_lc else 30)]
    logging.info("{0} datasets were de-trended".format(len(model.datasets)))
    if not len(datasets):
        raise RuntimeError("No datasets could be de-trended.")

    return model
