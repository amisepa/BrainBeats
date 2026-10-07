"""fractal_volatility - Fractal dimension of a 1-D signal by box counting.

Python port of BrainBeats functions/fractal_volatility.m (adapted from
fractalvol). The signal is rescaled into the unit square, boxes of width
2^-j counted down to the sampling resolution, and the dimension is the OLS
slope of log(count) vs log(1/width) after discarding scales whose local
slope deviates from the median by more than IQR/2.
Returns (dimension, standard_dev), dimension rounded to 3 decimals.

Copyright (C) Cedric Cannard, 2026 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import numpy as np

__all__ = ['fractal_volatility']


def fractal_volatility(data):
    """data: (n,) signal or (n, 2) [x y] matrix."""
    data = np.asarray(data, float)
    if data.ndim == 1:
        data = np.column_stack([np.arange(data.size, dtype=float), data])
    elif data.shape[0] < data.shape[1]:
        data = data.T                                   # MATLAB: rows >= cols
    if data.shape[1] == 1:
        data = np.column_stack([np.arange(data.shape[0], dtype=float), data[:, 0]])

    # normalize both axes to the unit square
    nrm = data.copy()
    for c in (0, 1):
        lo, hi = data[:, c].min(), data[:, c].max()
        nrm[:, c] = (data[:, c] - lo) / (hi - lo)

    # smallest box width stays above the sample spacing
    minwidth = np.ceil(np.abs(np.log2(np.diff(nrm[:, 0]).min())))
    minwidth = max(int(minwidth) - 1, 1)

    counts = np.zeros(minwidth)
    for j in range(1, minwidth + 1):
        width = 2.0 ** -j
        boxcount = 0
        xaxis_pos = 0.0
        while xaxis_pos < 1.0:
            indx = (nrm[:, 0] >= xaxis_pos) & (nrm[:, 0] < xaxis_pos + width)
            if abs(1.0 - xaxis_pos) == width:           # include the last sample
                indx[-1] = True
            col = nrm[indx, 1]
            if col.size == 1:
                boxcount += 1
            elif col.size > 1:
                raw = (col.max() - col.min()) / width + (col.min() % width)
                boxcount += int(np.ceil(raw))
            xaxis_pos += width
        counts[j - 1] = boxcount

    r = 2.0 ** -np.arange(1, minwidth + 1)
    log_r = np.log(r)
    log_n = np.log(counts)
    # local slopes; discard scales deviating from the median slope by > IQR/2
    s = -np.gradient(log_n, log_r)
    iqr = np.percentile(s, 75) - np.percentile(s, 25)
    keep = np.abs(s - np.median(s)) <= iqr / 2.0
    x2, y2 = log_r[keep], log_n[keep]

    X = np.column_stack([np.ones_like(x2), x2])
    beta = np.linalg.pinv(X) @ y2
    C = np.linalg.pinv(X.T @ X)
    e = y2 - X @ beta
    s2 = (e @ e) * C                            # MATLAB: e'*e * C (2x2)
    dimension = round(float(-beta[1]), 3)
    return dimension, float(np.sqrt(s2[1, 1]))