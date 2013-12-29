import os
import h5py
import numpy as np
from matplotlib import rc
rc("text", usetex=True)
import matplotlib.pyplot as pl

import kplr
client = kplr.API()

thin = 100
burnin = 10000

koi_ror = []
koi_ror_err = []
dfm_ror = []
dfm_ror_err = []
koi_b = []
koi_b_err = []
dfm_b = []
dfm_b_err = []
dfm_b_mean = []
dfm_b_map = []
bsamples = []
rorsamples = []

(kic_ids, koi_names, n_koi_tot, koi_num, koi_ror, koi_ror_err1, koi_ror_err2,
 koi_b, koi_b_err1, koi_b_err2) = zip(*[l.split()[:2]
                                        + map(int, l.split()[2:4])
                                        + map(float, l.split()[4:])
                                        for l in open("kois.txt", "r")])

for kic_id, n, i in zip(kic_ids, n_koi_tot, koi_num):
    fn = os.path.join("models", kic_id, "mcmc.h5")
    with h5py.File(fn, "r") as f:
        it = f.attrs["iteration"]
        ror = f["samples"][:, burnin:it:thin, 2+n*3+i].flatten()
        b = f["samples"][:, burnin:it:thin, 2+n*4+i].flatten()
        lp = f["lnprob"][:, burnin:it:thin].flatten()

    bsamples.append(b)
    rorsamples.append(ror)

    q = np.percentile(ror, (16, 50, 84))
    dfm_ror += [q[1]]
    dfm_ror_err += [(q[1]-q[0], q[2]-q[1])]

    q = np.percentile(b, (16, 50, 84))
    dfm_b += [q[1]]
    dfm_b_err += [(q[1]-q[0], q[2]-q[1])]
    dfm_b_mean.append(np.mean(b))
    dfm_b_map.append(b[np.argmax(lp)])

print("Running optimization")
bins = np.linspace(0, 1.2, 12)
counts = np.empty((len(bsamples), len(bins) - 1))
weights = np.empty((len(bsamples), len(bins) - 1))
for i, (b, ror) in enumerate(zip(bsamples, rorsamples)):
    counts[i, :], tmp = np.histogram(b, bins)
    weights[i, :], tmp = np.histogram(b, bins, weights=np.ones_like(ror))
weights /= np.mean(weights)


def lnlike(p):
    h = np.append(p, 1.0 - np.sum(p))
    if np.any(h < 0):
        return -np.inf
    vals = np.empty(len(bsamples))
    for i, b in enumerate(bsamples):
        vals[i] = np.log(np.sum(weights[i]*h)) - np.sum(counts[i])
    ll = np.sum(vals)
    return ll


p0 = np.ones(len(bins) - 1)
p0 /= np.sum(p0)
p0 = p0[:-1]

import emcee
ndim, nwalkers = len(p0), 50
pos = [p0 + 1e-8*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike)
pos, lp, state = sampler.run_mcmc(pos, 1000)
sampler.reset()
sampler.run_mcmc(pos, 1000)

for i in range(sampler.chain.shape[-1]):
    pl.plot(sampler.chain[:, :, i].T)
    pl.savefig("time/{0}.png".format(i))
    pl.clf()

samples = np.hstack((sampler.flatchain,
                     np.atleast_2d(1-np.sum(sampler.flatchain, axis=1)).T))

b = np.array(zip(bins[:-1], bins[1:])).T
for ind in np.random.randint(len(samples), size=100):
    v = samples[ind]/(bins[1] - bins[0])
    pl.plot(b, np.array(zip(v, v)).T, "k", alpha=0.1)

v = np.median(samples, axis=0)/(bins[1] - bins[0])
pl.plot(b.T.flatten(), np.array(zip(v, v)).flatten(), "k", lw=2)

pl.hist(koi_b, bins, histtype="step", normed=True, label="koi", lw=2)
pl.hist(dfm_b, bins, histtype="step", normed=True, label="med", lw=2)
pl.hist(dfm_b_mean, bins, histtype="step", normed=True, label="mean", lw=2)
pl.hist(dfm_b_map, bins, histtype="step", normed=True, label="map", lw=2)
pl.legend()

pl.savefig("sampling.png")

assert 0
pl.figure(figsize=(10, 10))
pl.errorbar(koi_ror, dfm_ror, xerr=(np.abs(koi_ror_err2),
                                    koi_ror_err1),
            yerr=np.array(dfm_ror_err).T, fmt=".k")
pl.plot([0, 1], [0, 1], ":r")

for i, kic_id in enumerate(kic_ids):
    pl.gca().annotate(kic_id, xy=(koi_ror[i], dfm_ror[i]), xycoords="data",
                      ha="left", va="top")

pl.xlabel(r"$(r/R)_\mathrm{KOI}$")
pl.ylabel(r"$(r/R)_\mathrm{DFM}$")
pl.xlim(0, 1)
pl.ylim(0, 1)
pl.savefig("ror.png")

pl.xlim(0, 0.1)
pl.ylim(0, 0.1)
pl.savefig("ror_zoom.png")

pl.clf()
pl.errorbar(koi_b, dfm_b, xerr=(np.abs(koi_b_err2), koi_b_err1),
            yerr=np.array(dfm_b_err).T, fmt=".k")
pl.plot([0, 1], [0, 1], ":r")
pl.xlabel(r"$b_\mathrm{KOI}$")
pl.ylabel(r"$b_\mathrm{DFM}$")
pl.xlim(0, 2)
pl.ylim(0, 2)
pl.savefig("b.png")
