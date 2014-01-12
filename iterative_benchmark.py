#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = []

import time
import kplr
import numpy as np

from kois.kois import Model
from kois.lightcurve.ld import QuadraticLimbDarkening

P = 150.0
t = np.concatenate([i*P + np.arange(0.0, 10.0, 0.5/24.) for i in range(100)])

ldp = QuadraticLimbDarkening(0.3, 0.1)
model = Model("test", ldp, tol=1e-20, max_depth=10000)
model.add_koi(P, 5.0, 13.0 / 24.0, 0.5, 0.8)

print(model.durations)
texp = kplr.EXPOSURE_TIMES[1]/86400.0

for tau in [0.05, 0.1, 0.5, 1.0]:
    model.durations[0] = tau
    strt = time.time()
    lc = model.get_light_curve(t, texp=texp)
    print(tau, time.time() - strt)
