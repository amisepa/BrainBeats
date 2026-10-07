"""estimate_pat - Pulse arrival time: delay between each ECG R-peak and the
following PPG pulse fiducial.

Python port of BrainBeats functions/estimate_pat.m (v1.6). Each PPG beat is
paired with the closest PRECEDING R-peak within patRange ms (default
[50 600]); the median PAT shifts PPG beats back to the heartbeats.

Copyright (C) Cedric Cannard, 2026 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import numpy as np

__all__ = ["estimate_pat"]


def estimate_pat(rpeaks, ppgpeaks, fs, pat_range=(50.0, 600.0)):
    rpeaks = np.sort(np.asarray(rpeaks, float).ravel())
    ppg = np.asarray(ppgpeaks, float).ravel()
    pat = np.full(ppg.size, np.nan)
    iR = np.searchsorted(rpeaks, ppg, side="left") - 1   # closest preceding
    ok = (iR >= 0)
    d = np.full(ppg.size, np.nan)
    d[ok] = (ppg[ok] - rpeaks[iR[ok]]) / fs * 1000.0
    inr = (d >= pat_range[0]) & (d <= pat_range[1])
    pat[inr] = d[inr]
    n_paired = int(np.isfinite(pat).sum())
    info = {"median": float(np.nanmedian(pat)) if n_paired else float("nan"),
            "iqr": (float(np.nanpercentile(pat, 75) - np.nanpercentile(pat, 25))
                    if n_paired else float("nan")),
            "sd": float(np.nanstd(pat)) if n_paired else float("nan"),
            "nPaired": n_paired, "nPPG": int(ppg.size)}
    if n_paired < 0.5 * info["nPPG"]:
        import warnings
        warnings.warn(f"estimate_pat: only {n_paired} of {info['nPPG']} PPG beats "
                      "could be paired with an R-peak: check that both signals "
                      "are aligned.")
    return info["median"], pat, info
