#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging
import cPickle as pickle

import kois

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build a KOI model")
    parser.add_argument("koi", type=int, help="The KOI ID")
    parser.add_argument("-o", "--outdir", default=None,
                        help="The directory where the output files will be "
                        "written")
    parser.add_argument("--plots", action="store_true",
                        help="Make some basic figures")
    args = parser.parse_args()

    # Set up the output directory.
    outdir = args.outdir
    if outdir is None:
        outdir = "models"
    outdir = os.path.join(outdir, "koi-{0}".format(args.koi))
    try:
        os.makedirs(outdir)
    except os.error:
        logging.info("Output directory '{0}' exists".format(outdir))

    model = kois.load_system(args.koi)
    pickle.dump(model, open(os.path.join(outdir, "model.pkl"), "wb"), -1)

    if args.plots:
        import numpy as np
        import matplotlib.pyplot as pl

        for i, (P, t0) in enumerate(zip(model.periods, model.epochs)):
            hp = 0.5 * P
            pl.clf()
            [pl.plot((d.time-t0+hp) % P - hp, d.flux, ".k", alpha=0.5)
             for d in model.datasets]

            t = np.concatenate([d.time for d in model.datasets])
            pl.plot((t-t0+hp) % P - hp, model.get_light_curve(t), ".r", ms=3)

            pl.xlim(-1, 1)
            pl.title("KOI {0}.{1:02d}, period = {2} days"
                     .format(args.koi, i+1, P))
            pl.savefig(os.path.join(outdir,
                                    "light-curve-{0}.{1:02d}.png"
                                    .format(args.koi, i+1)))

# nwalkers, ndim = 64, len(model.vector)
# p0 = [model.vector + 1e-6 * np.random.randn(ndim) for i in range(nwalkers)]
# sampler = emcee.EnsembleSampler(nwalkers, ndim, model)
# for pos, lp, state in sampler.sample(p0, iterations=500):
#     print(lp)
