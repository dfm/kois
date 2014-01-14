#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = []

import time
import kplr
import numpy as np
import matplotlib.pyplot as pl

from kois.kois import Model
from kois.lightcurve.ld import QuadraticLimbDarkening

P = 150.0
# t = np.concatenate([i*P + np.arange(0.0, 10.0, 0.5/24.) for i in range(100)])

ldp = QuadraticLimbDarkening(0.3, 0.1)
model = Model("test", ldp, tol=1e-12, max_depth=50)
# model.add_koi(P, 0.0, 12.0 / 24.0, 0.04, 0.7)
model.add_koi(P, 0.0, 0.4 / 24.0, 0.04, 0.7)

texp = kplr.EXPOSURE_TIMES[1]/86400.0
print(texp)

t = np.linspace(-1, 1, 10013)
# t = 0.02 - 0.04 * np.sort(np.random.rand(200))

lcs = []
for tol, md, fmt in [(0.0, 1, "-r"),
                     (0.0, 2, "--k"),
                     (0.0, 3, "--k"),
                     ]:
    print(fmt)
    model.tol = tol
    model.max_depth = md

    strt = time.time()
    lcs.append(model.get_light_curve(t, texp=texp))
    print(time.time() - strt)

pl.plot(t, lcs[2] - lcs[0], "k")
pl.plot(t, lcs[2] - lcs[1], "r")

pl.xlim(-0.03, 0.03)
# pl.ylim(0.9987, 1.0002)
pl.savefig("diff.png")



# print(model.durations)

# for tau in [0.05, 0.1, 0.5, 1.0]:
#     model.durations[0] = tau
#     strt = time.time()
#     lc = model.get_light_curve(t, texp=texp)
#     print(tau, time.time() - strt)
