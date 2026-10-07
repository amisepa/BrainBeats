"""baseline_regression - regression-based baseline correction of epoched data.

Port of functions/baseline_regression.m (Alday, 2019). Instead of subtracting
each trial's baseline, the baseline mean is used as a trial-level regressor;
the slope is estimated at each channel x time point and each trial is corrected
by that slope times its centered baseline. The trial average is unchanged,
which a plain subtraction is not for HEPs (the pre-R window holds the previous
cardiac cycle; Steinfath et al. 2026 recommend no baseline correction, or this
regression form).

    data, beta, bl = baseline_regression(data, times, bl_win)
    data, beta, bl = baseline_regression(..., groups)

data : (channels x times x trials) epochs
times: epoch time vector in ms
bl_win: baseline window [start end] in ms
groups: optional per-trial condition labels -> within-group slopes (ANCOVA)

Copyright (C) Cedric Cannard, 2026 -- Python port, 2026-10.
"""
from __future__ import annotations

import numpy as np

__all__ = ['baseline_regression']


def baseline_regression(data, times, bl_win, groups=None):
    data = np.asarray(data, dtype=float)
    if data.ndim == 2:                      # (times x trials) single channel
        data = data[None]
    n_chan, n_times, n_trials = data.shape
    times = np.asarray(times, dtype=float).ravel()
    bl_idx = (times >= bl_win[0]) & (times <= bl_win[1])
    if not bl_idx.any():
        raise ValueError(f'baseline_regression: no time point in the baseline '
                         f'window [{bl_win[0]} {bl_win[1]}] ms.')
    if n_trials < 3:
        raise ValueError('baseline_regression: at least 3 trials are needed.')

    if groups is None:
        g = np.ones(n_trials, dtype=int)
        groups_list = [1]
    else:
        labels = list(dict.fromkeys(np.asarray(groups).ravel().tolist()))
        lut = {lab: i + 1 for i, lab in enumerate(labels)}
        g = np.array([lut[v] for v in np.asarray(groups).ravel().tolist()])
        groups_list = list(range(1, len(labels) + 1))
    n_g = len(groups_list)

    bl = data[:, bl_idx, :].mean(axis=1)          # channels x trials
    blc = bl - bl.mean(axis=1, keepdims=True)     # centered over ALL trials

    beta = np.full((n_chan, n_times), np.nan)
    # X is the same for every channel: intercepts per group + centered baseline
    G = np.zeros((n_trials, n_g))
    for j, gi in enumerate(groups_list):
        G[g == gi, j] = 1.0
    # pinv per channel only through the baseline column (channels differ)
    for i in range(n_chan):
        X = np.column_stack([G, blc[i]])
        Y = data[i].T                             # trials x times
        coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
        beta[i] = coef[-1]

    corrected = data - beta[:, :, None] * blc[:, None, :]
    return corrected, beta, bl