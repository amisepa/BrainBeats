"""compute_fe - Fuzzy entropy (FuzzyEn) of a univariate signal.

Python port of BrainBeats functions/compute_fe.m. The signal is z-scored
(NaNs ignored), so `r` is a fraction of its SD. Similarity between embedding
vectors is exp(-d^n / r), with d the Chebyshev distance (self-matches
excluded). Returns (fe, p) with p the global similarity in dimensions m and
m+1; fe = ln(p[0]/p[1]) rounded to 3 decimals (MATLAB parity).

Reference: Azami & Escudero 2016; Chen et al. 2007.

Copyright (C) Cedric Cannard, 2026 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import numpy as np

__all__ = ['compute_fe']


def _similarity_counts(xm, m, r, n):
    """Global similarity p for dimensions m and m+1 (compute_fe.m:55-66)."""
    N = xm.size
    n_emb = N - m                     # embedding vectors of dimension m
    # lagged copies, one row per tap (m+1 rows), MATLAB xMat(i,:) =
    # signal(i : N-m+i-1) -> python xMat[i-1, :] = signal[i-1 : N-m+i-1]
    xMat = np.empty((m + 1, n_emb))
    for i in range(m + 1):
        xMat[i] = xm[i:i + n_emb]
    p = np.zeros(2)
    for k_off, k in enumerate((m, m + 1)):       # k = embedding size in rows
        count = np.zeros(n_emb)
        tmp = xMat[:k]                           # xMat(1:k, :)
        # forward pairs only (MATLAB i = 1:N-k excludes self-matches and
        # double counting): dist(i, j) = Chebyshev(vec i, vec j), j > i
        dists = np.max(np.abs(tmp[:, :, None] - tmp[:, None, :]), axis=0)
        # dists[i, j] = max over taps |tmp[:, i] - tmp[:, j]|
        for i in range(N - k):
            d = dists[i, i + 1:n_emb]
            df = np.exp(-(d ** n) / r)
            count[i] = df.sum() / n_emb
        p[k_off] = count.sum() / n_emb
    return p


def compute_fe(signal, m=2, r=0.15, n=2, tau=1, useGPU=False):
    """Fuzzy entropy of a 1-D signal (row or column), MATLAB-compatible.

    m: embedding dimension (2); r: tolerance as a fraction of the signal's
    SD (0.15); n: fuzzy power (2); tau: time lag, > 1 decimates by tau (1).
    Returns (fe, p): fe = ln(p(m)/p(m+1)) rounded to 3 decimals (NaN-safe),
    p = [p_m, p_m+1].
    """
    x = np.asarray(signal, float).ravel()
    if tau > 1:
        x = x[::tau]                       # MATLAB downsample: every tau-th
    if np.all(np.isnan(x)):
        return float('nan'), np.array([np.nan, np.nan])
    mu = np.nanmean(x)
    sd = np.nanstd(x)
    x = (x - mu) / sd                      # z-score (NaN-omit)
    if not np.isfinite(sd) or sd == 0:
        return float('nan'), np.array([np.nan, np.nan])
    xm = np.nan_to_num(x)                  # MATLAB: NaN stays in place and
    # propagates through distances as NaN -> 'omitnan' in sum treats them as
    # missing; after a clean z-score NaNs only exist at edge artefacts. The
    # MATLAB code sums with 'omitnan' over NaN distances: a NaN in ANY tap
    # makes that pair uncountable. Mirror it by setting such pairs to 0
    # similarity is NOT equal -- replicate: distances involving NaN taps give
    # NaN df -> sum(df,'omitnan') drops them from the SUM but keeps the
    # denominator N-m. Emulate with an output-ignoring nan-sum.
    p = _similarity_counts_nanaware(xm, m, r, n, isnan=np.isnan(x))
    with np.errstate(divide='ignore', invalid='ignore'):
        fe = np.log(p[0] / p[1])
    fe = round(float(fe), 3)
    return fe, p


def _similarity_counts_nanaware(xm, m, r, n, isnan):
    """_similarity_counts with MATLAB's sum(...,'omitnan') behaviour."""
    N = xm.size
    n_emb = N - m
    xMat = np.empty((m + 1, n_emb))
    for i in range(m + 1):
        xMat[i] = xm[i:i + n_emb]
    nanMat = np.empty((m + 1, n_emb), bool)
    for i in range(m + 1):
        nanMat[i] = isnan[i:i + n_emb]
    p = np.zeros(2)
    for k_off, k in enumerate((m, m + 1)):
        count = np.zeros(n_emb)
        tmp = xMat[:k]
        tmpnan = nanMat[:k]
        # any NaN in either embedding vector -> that pair contributes NaN
        bad = tmpnan.any(axis=0)
        dists = np.max(np.abs(tmp[:, :, None] - tmp[:, None, :]), axis=0)
        with np.errstate(invalid='ignore'):
            D = np.where(bad[:, None] | bad[None, :], np.nan, dists)
        for i in range(N - k):
            d = D[i, i + 1:n_emb]
            with np.errstate(invalid='ignore'):
                df = np.exp(-(d ** n) / r)
            count[i] = np.nansum(df) / n_emb        # omitnan
        p[k_off] = count.sum() / n_emb
    return p