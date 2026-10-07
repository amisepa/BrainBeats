"""resample_NN - Resample an NN interval series on a regular time grid.

Python port of BrainBeats functions/resample_NN.m. The grid runs from
NN_times[0] to NN_times[-1] at `sf` Hz; interpolation 'cub'/'spline' =
cubic spline, 'lin'/'linear' = linear, 'pchip' = shape-preserving cubic.

Copyright (C) Cedric Cannard, 2024 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator, interp1d

__all__ = ['resample_NN']


def resample_NN(NN_times, NN, sf, interp_method):
    """Returns (NN_resamp, t_resamp) on the regular grid between the ends."""
    t = np.asarray(NN_times, float).ravel()
    x = np.asarray(NN, float).ravel()
    if t.size != x.size:
        raise ValueError('resample_NN: NN and NN_times must have the same length')
    t_resamp = np.arange(t[0], t[-1] + 1.0 / sf, 1.0 / sf)
    t_resamp = t_resamp[t_resamp <= t[-1]]          # MATLAB colon: end-inclusive
    m = interp_method.lower()
    if m in ('cub', 'spline'):
        f = CubicSpline(t, x)
    elif m in ('lin', 'linear'):
        f = interp1d(t, x, kind='linear')
    elif m == 'pchip':
        f = PchipInterpolator(t, x)
    else:
        raise ValueError(f"resample_NN: unknown interp_method '{interp_method}'. "
                         "Use 'cub'/'spline', 'lin'/'linear' or 'pchip'.")
    return f(t_resamp), t_resamp