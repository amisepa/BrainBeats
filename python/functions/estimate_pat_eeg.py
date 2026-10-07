"""estimate_pat_eeg - Pulse arrival time (PAT) of a PPG estimated from the
cardiac field artifact (CFA) of the EEG, when no ECG was recorded.

Python port of BrainBeats functions/estimate_pat_eeg.m (v1.6). The QRS field
reaches the scalp with no delay while the PPG pulse arrives ~200-450 ms
later: the GFP peak of the median EEG average around the PPG beats between
-650 and -100 ms gives the PAT. Kept only if peak/median GFP >= 2, else NaN.

Copyright (C) Cedric Cannard, 2026 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy.signal import butter, filtfilt

__all__ = ["estimate_pat_eeg"]


def estimate_pat_eeg(X, fs, ppgbeats, search_win=(-650.0, -100.0), min_ratio=2.0):
    X = np.asarray(X, float) - np.asarray(X, float).mean(axis=0, keepdims=True)
    b, a = butter(2, np.asarray([5.0, 30.0]) / (fs / 2.0), btype="band")
    X = filtfilt(b, a, X, axis=1)
    lags = np.arange(int(round(-0.7 * fs)), int(round(0.2 * fs)) + 1)
    t = lags / fs * 1000.0
    ppg = np.round(np.asarray(ppgbeats, float).ravel()).astype(np.int64)
    ppg = ppg[(ppg + lags[0] >= 1) & (ppg + lags[-1] <= X.shape[1])]
    n_b = ppg.size
    info = {"pat": float("nan"), "ratio": float("nan"), "nBeats": int(n_b),
            "gfp": np.empty(0), "times": t}
    if n_b < 30:
        return float("nan"), info
    # median across beats (robust to artifacts), then GFP
    pos = (lags[:, None] + ppg[None, :]) - 1                  # 0-based: lags x beats
    E = X[:, pos]                                             # ch x lags x beats
    gfp = np.median(E, axis=2).std(axis=0, ddof=1)            # GFP over channels
    w = (t >= search_win[0]) & (t <= search_win[1])
    i_max = int(np.argmax(gfp[w]))
    info["gfp"] = gfp
    info["ratio"] = float(gfp[w][i_max] / np.median(gfp))
    if info["ratio"] >= min_ratio:
        info["pat"] = float(-t[w][i_max])
    else:
        warnings.warn("estimate_pat_eeg: cardiac field artifact peak not clear "
                      f"(ratio {info['ratio']:.1f} < {min_ratio:g}).")
    return info["pat"], info