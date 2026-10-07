"""compute_asymmetry - Alpha asymmetry on symmetric electrode pairs:
ln(alpha_left + eps) - ln(alpha_right + eps).

Python port of BrainBeats functions/compute_asymmetry.m (Allen et al. 2004;
Smith et al. 2017). Each LEFT electrode (Y > tol) is paired with the RIGHT
electrode (Y < -tol) closest to its mirror position XYZ .* [1 -1 1] (within
10% of the head radius); midline ('z') labels are excluded. Returns
(asy, pairLabels, pairNums) with pairNums [left right] 0-based indices.

Copyright (C) Cedric Cannard, 2023 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import numpy as np

__all__ = ['compute_asymmetry']


def _xyz_from_chanlocs(chanlocs):
    """EEGLAB chanlocs (dict-like rows) -> (nChan, 3) X/Y/Z array."""
    cols = {'X': [], 'Y': [], 'Z': []}
    for c in np.asarray(chanlocs, dtype=object).ravel():
        d = c if isinstance(c, dict) else None
        if d is None:
            cols['X'].append(np.nan); cols['Y'].append(np.nan); cols['Z'].append(np.nan)
            continue
        for key in ('X', 'Y', 'Z'):
            v = d.get(key if key in d else key.lower(), np.nan)
            try:
                cols[key].append(float(v))
            except (TypeError, ValueError):
                cols[key].append(np.nan)
    return np.column_stack([cols['X'], cols['Y'], cols['Z']])


def _labels_from_chanlocs(chanlocs):
    return [str(c['labels'] if isinstance(c, dict) else getattr(c, 'labels', c))
            for c in np.asarray(chanlocs, dtype=object).ravel()]


def compute_asymmetry(alpha_pwr, norm, chanlocs, vis=False, tot_pwr=None):
    """alpha_pwr: mean alpha PSD per channel in uV^2/Hz (NOT dB).
    norm: divide by each channel's total power (tot_pwr) first.
    Returns (asy (nPairs,), pairLabels (nPairs,), pairNums (nPairs, 2))."""
    if norm and (tot_pwr is None or (np.asarray(tot_pwr).size == 0)):
        raise ValueError("compute_asymmetry: 'tot_pwr' is required to normalize "
                         "asymmetry.")
    labels = _labels_from_chanlocs(chanlocs)
    nChan = len(labels)
    pairNums = np.full((nChan, 2), np.nan)
    pairLabels = [None] * nChan
    XYZ = _xyz_from_chanlocs(chanlocs)
    if np.isnan(XYZ).any(axis=None) and np.isnan(XYZ).all():
        print('compute_asymmetry warning: no electrode positions; no pairs.')
        return np.array([]), [], np.empty((0, 2))
    r = float(np.nanmedian(np.sqrt((XYZ ** 2).sum(axis=1))))
    tol = 0.02 * r
    for iChan in np.flatnonzero(XYZ[:, 1] > tol):
        mirrored = XYZ[iChan] * np.array([1.0, -1.0, 1.0])
        d = np.sqrt(((XYZ - mirrored) ** 2).sum(axis=1))
        d[iChan] = np.inf
        match = int(np.argmin(d))
        dmin = d[match]
        if dmin < 0.1 * r and XYZ[match, 1] < -tol:
            pairNums[iChan] = [iChan, match]
            pairLabels[iChan] = f'{labels[iChan]} {labels[match]}'

    ok = [i for i in range(nChan) if pairLabels[i] is not None]
    pairNums = pairNums[ok] if ok else np.empty((0, 2))
    pairLabels = [pairLabels[i] for i in ok]

    # drop pairs touching midline ('z') electrodes, if any slipped through
    keep = [i for i, lab in enumerate(pairLabels) if 'z' not in lab.lower()]
    pairNums = pairNums[keep] if keep else np.empty((0, 2))
    pairLabels = [pairLabels[i] for i in keep]
    if len(pairNums) != len(pairLabels):
        import warnings
        warnings.warn('compute_asymmetry: pair count mismatch (labels vs indices).')

    alpha_pwr = np.asarray(alpha_pwr, float).ravel()
    if norm:
        alpha_pwr = alpha_pwr / np.asarray(tot_pwr, float).ravel()
    alpha_log = np.log(alpha_pwr + np.finfo(float).eps)

    asy = np.full(len(pairLabels), np.nan)
    for iPair in range(len(pairLabels)):
        left, right = pairNums[iPair]
        asy[iPair] = alpha_log[int(left)] - alpha_log[int(right)]
    return asy, pairLabels, pairNums.astype(int)