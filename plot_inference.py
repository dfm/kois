import numpy as np
import matplotlib.pyplot as pl

fn = "rdist_samples.txt"
bins = []
for line in open(fn):
    bins += line.strip().strip("]").split("[")[-1].split()
    if "]" in line:
        break

bins = np.array(bins, dtype=float)
samples = np.loadtxt(fn, skiprows=10007)[:, 1:-1]

samples = np.concatenate((samples, 1-np.atleast_2d(np.sum(samples, axis=1)).T),
                         axis=1)

b = np.array(zip(bins[:-1], bins[1:])).T
v = np.median(samples, axis=0)
pl.plot(b.T.flatten(), np.array(zip(v, v)).flatten())
pl.savefig("ror_inference.png")
