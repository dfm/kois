#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from setuptools import setup, Extension

includes = [
    "kois",
    numpy.get_include(),
]

kois = Extension("kois._kois", ["kois/_kois.c", "kois/lightcurve.c",
                                "kois/quad.cpp"],
                 include_dirs=includes)

setup(
    name="kois",
    url="https://github.com/dfm/kois",
    version="0.0.0",
    author="Dan Foreman-Mackey",
    author_email="danfm@nyu.edu",
    description="",
    long_description="",
    packages=["kois"],
    ext_modules=[kois],
    install_requires=[
        "numpy",
    ],
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
