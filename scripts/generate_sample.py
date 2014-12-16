#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

__all__ = ["get_stars", "get_kois"]

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

d = os.path.dirname
data_dir = os.path.join(d(d(os.path.abspath(__file__))), "raw_data")
output_dir = os.path.join(d(d(os.path.abspath(__file__))), "data")


def get_stars():
    # Load the revised KIC properties.
    fn = os.path.join(data_dir, "huber-kic-join-2014-12-16.tsv.gz")
    kic = pd.read_csv(fn, sep="\t", compression="gzip")

    # Make the magnitude cut.
    m = kic.kic_kepmag < 13.5

    # Make the T_eff cut.
    m &= (4100 < kic.Teff) & (kic.Teff < 6100)

    # Make the log_g cut.
    m &= (4 < kic["log(g)"]) & (kic["log(g)"] < 4.9)

    return kic[m]


def get_kois(stars, min_period=100.0):
    kois = pd.read_csv(os.path.join(data_dir, "kois-2014-12-16.tsv"),
                       sep="\t", skiprows=155)

    # Select only candidates.
    m = kois.koi_disposition == "CANDIDATE"
    m |= kois.koi_disposition == "CONFIRMED"

    # Select only long-period candidates.
    m &= kois.koi_period > min_period

    # Select only KOIs in the stellar sample.
    joined = pd.merge(kois[m], stars, left_on="kepid", right_on="KIC")

    return joined


if __name__ == "__main__":
    stars = get_stars()
    kois = get_kois(stars)

    # Save the output files.
    try:
        os.makedirs(output_dir)
    except os.error:
        pass
    kois.to_hdf(os.path.join(output_dir, "kois.h5"), "kois")
    stars.to_hdf(os.path.join(output_dir, "stars.h5"), "stars",
                 complib="zlib", complevel=9)
