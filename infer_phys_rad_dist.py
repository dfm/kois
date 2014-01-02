#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = []

import os
import h5py
import emcee
import sqlite3
import numpy as np
import multiprocessing
from scipy.misc import logsumexp
from scipy.linalg import cho_factor, cho_solve
from collections import defaultdict


def load_dressing_stars():
    fn = os.path.join("kois", "resources", "dressing_stars.txt")
    lines = open(fn, "r").readlines()
    catalog = {}
    for l in lines[28:]:
        cols = l.split()
        catalog[cols[0]] = np.array(cols[1:], dtype=float)
    return catalog


def load_cool_kois():
    with sqlite3.connect("results/kois.db") as conn:
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        c.execute("select * from kois where kplr_ror is not NULL")
        return c.fetchall()


def load_koi_samples(row, size=500):
    fn = os.path.join("results", row[str("kepoi_name")], "mcmc.h5")
    with h5py.File(fn, "r") as f:
        i = f.attrs["iteration"]
        chain = f["samples"][:, row[str("burnin")]:i, :]
    w = np.random.randint(chain.shape[0], size=size)
    i = np.random.randint(chain.shape[1], size=size)
    return chain[w, i, 6]


print("Loading catalog")
dressing_catalog = load_dressing_stars()
kois = load_cool_kois()
nz = 100
full_catalog = defaultdict(lambda: [None, np.empty((0, nz))])
for koi in kois:
    # Draw some samples from the Dressing stellar radius.
    k = koi[str("kepid")]
    if full_catalog[k][0] is None:
        row = dressing_catalog[str(k)]
        R = row[3:6]
        Rs = R[0] + 0.5*sum(R[1:])*np.random.randn(200)
        full_catalog[k][0] = Rs

    # Grab some samples from the KOI MCMC chain.
    full_catalog[k][1] = np.concatenate(
        (full_catalog[k][1],
            np.atleast_2d(load_koi_samples(koi, nz))))

# Set up the histogram.
bins = np.linspace(0, 1.0, 40)
values = np.random.rand(len(bins) - 1)
values /= np.sum(values)
values = values[:-1]

# Pre-compute the GP kernel and factorize.
a2 = 1e-4
l2 = 0.01 ** 2
x = 0.5*(bins[:-1]+bins[1:])
K = a2*np.exp(-(x[:, None] - x[None, :]) ** 2 / l2)
factor = cho_factor(K)


def lnlike(v):
    ll = 0.0
    for R, z in full_catalog.itervalues():
        r = R[:, None, None] * z[None, :, :]
        inds = np.digitize(r.flatten(), bins)
        pr = v[inds].reshape(r.shape)*(z*z)[None, :, :]
        pr = np.sum(pr, axis=2)
        if np.any(pr == 0):
            return -np.inf
        pr = np.sum(np.log(pr), axis=1)
        if not np.all(np.isfinite(pr)):
            return -np.inf
        ll += logsumexp(pr)
    return ll


def lnprior(v):
    if np.any(v < 0.0):
        return -np.inf
    return 0.0
    # return -0.5 * np.dot(v, cho_solve(factor, v))


def lnprob(p):
    v = np.append(p, 1.0 - np.sum(p))
    lp = lnprior(v)
    if not np.isfinite(lp):
        return -np.inf
    ll = lnlike(v)
    if not np.isfinite(ll):
        return -np.inf
    return ll + lp


print("Computing initial ln-prob")
print(lnprob(values))
ndim, nwalkers = len(values), 100
p0 = [np.abs(values+1e-8*np.random.randn(ndim)) for i in range(nwalkers)]
pool = multiprocessing.Pool()
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)

outfn = "rdist_samples.txt"
open(outfn, "w").close()
for i, (p, lp, s) in enumerate(sampler.sample(p0, iterations=5000)):
    if i % 10 == 0:
        print(i, np.mean(sampler.acceptance_fraction))
        with open(outfn, "a") as f:
            for l, v in zip(p, lp):
                f.write("{0} ".format(i) + " ".join(map("{0}".format, l))
                        + " {0}\n".format(v))
