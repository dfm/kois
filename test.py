#!/usr/bin/env python
# -*- coding: utf-8 -*-

import emcee
import numpy as np
import matplotlib.pyplot as pl

import kois

model = kois.load_system("117")

nwalkers, ndim = 64, len(model.vector)
p0 = [model.vector + 1e-6 * np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, model)
for pos, lp, state in sampler.sample(p0, iterations=500):
    print(lp)


# for i, (P, t0) in enumerate(zip(model.periods, model.epochs)):
#     hp = 0.5 * P
#     pl.clf()
#     [pl.plot((d.time-t0+hp) % P - hp, d.flux, ".") for d in model.datasets]
#     t = np.linspace(-1, 1, 1000)
#     pl.plot(t, model.get_light_curve(t+t0), "k")
#     pl.xlim(-1, 1)
#     pl.savefig("data.{0}.png".format(i))
