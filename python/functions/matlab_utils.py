"""MATLAB-semantics helpers: 1-based indices, MATLAB rounding, FDR, Grubbs."""
from __future__ import annotations

import numpy as np

__all__ = ['matlab_round', 'fdr_bh', 'grubbs_outliers', 'percentile_matlab']


def matlab_round(x):
    """MATLAB round: half away from zero (numpy rounds half to even)."""
    x = np.asarray(x, dtype=float)
    return np.where(x >= 0, np.floor(x + 0.5), np.ceil(x - 0.5))


def fdr_bh(p):
    """Benjamini-Hochberg FDR-adjusted p-values (MATLAB fdr_bh port).

    Same size as p; NaN ignored (stays NaN). Vectorised form of
    compute_hep_tf.m:284-292.
    """
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    good = np.isfinite(p)
    if not good.any():
        return q
    pv = p[good]
    order = np.argsort(pv, kind='stable')
    ps = pv[order]
    m = ps.size
    # cummin from the largest p backwards (MATLAB flipud(cummin(flipud(...))))
    rev = ps[::-1] * m / np.arange(m, 0, -1)
    adj = np.minimum.accumulate(rev)[::-1]
    adj = np.minimum(adj, 1.0)
    q.flat[np.flatnonzero(good)[order]] = adj
    return q


def grubbs_outliers(x, alpha: float = 0.05):
    """MATLAB isoutlier(x, 'grubbs'): two-sided Grubbs test, iterated.

    G = |x - mean| / s (s = sample sd, N-1) against
    Gcrit = (n-1)/sqrt(n) * sqrt(tc^2 / (n - 2 + tc^2)) with
    tc = t-inv at 1 - alpha/(2n) (NIST two-sided Grubbs, alpha SPLIT OVER THE
    TWO TAILS AND n CANDIDATES; using plain alpha/2 here would make Gcrit
    ~1.96 for n ~ 200, i.e. EVERY normal point an outlier -- that was caught
    by test_grubbs_matches_isoutlier_semantics). MATLAB isoutlier's 'grubbs'
    method documents the same test at the requested significance level.
    Iterate: remove the detected point, recompute, stop when none exceeds.
    Returns a boolean mask, True = outlier.
    """
    x = np.asarray(x, dtype=float)
    mask = np.zeros(x.shape, dtype=bool)
    good = np.isfinite(x)
    v = x[good]
    if v.size < 3:
        return mask
    idx_good = np.flatnonzero(good)
    while v.size >= 3:
        mu = v.mean()
        s = v.std(ddof=1)
        if s == 0:
            break
        i = int(np.argmax(np.abs(v - mu)))
        tmax = abs(v[i] - mu) / s
        n = v.size
        tcrit = _tcrit_two_sided(alpha, n)
        crit = ((n - 1) / np.sqrt(n)) * np.sqrt(tcrit ** 2 / (n - 2 + tcrit ** 2))
        if tmax > crit:
            mask[idx_good[i]] = True
            # remove from the working set
            keep = np.ones(v.size, dtype=bool)
            keep[i] = False
            v = v[keep]
            idx_good = idx_good[keep]
        else:
            break
    return mask


def _tcrit_two_sided(alpha: float, n: int) -> float:
    """Two-sided t critical value for Grubbs: t(1 - alpha/(2n), n-2)."""
    from scipy.stats import t as tdist
    return float(tdist.ppf(1 - alpha / (2.0 * n), n - 2))


def percentile_matlab(x, p):
    """MATLAB prctile: linear interpolation on (100*(i-0.5))/n plotting positions.

    numpy's default 'linear' method uses (i-1)/(n-1); MATLAB differs at the
    extremes. Used by the adaptive HEP window.
    """
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    xs = np.sort(x)
    n = xs.size
    pos = (np.asarray(p, float) / 100.0) * n + 0.5
    pos = np.clip(pos, 1, n)
    lo = np.floor(pos).astype(int) - 1
    hi = np.minimum(lo + 1, n - 1)
    frac = pos - np.floor(pos) if False else (pos - (lo + 1))
    # pos is 1-based; interpolate between xs[lo] and xs[hi]
    out = xs[lo] + (pos - (lo + 1)) * (xs[hi] - xs[lo])
    return float(out) if np.ndim(p) == 0 else out