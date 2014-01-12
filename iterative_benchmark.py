#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = []

import numpy as np

P = 150.0
t = np.concatenate([i*P + np.arange(0.0, 10.0, 0.5/24.) for i in range(4)])
print(t)
