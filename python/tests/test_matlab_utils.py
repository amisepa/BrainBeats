"""Unit tests for the MATLAB-semantics primitives (no MATLAB needed)."""
import sys
import numpy as np
import pytest

sys.path.insert(0, r'C:\Users\ccann\Documents\MATLAB\BrainBeats\python')
from functions.matlab_utils import (matlab_round, fdr_bh, grubbs_outliers,
                                    percentile_matlab)  # noqa: E402
from functions.compute_hep_tf import (valid_beats, _tf_stats, morlet_kernels_bb,
                                      surrogate_trains, surrogate_stats,
                                      compute_hep_tf)  # noqa: E402


def test_matlab_round():
    assert np.allclose(matlab_round(np.array([-2.5, 2.5, 3.5, -3.5])),
                       [-3, 3, 4, -4])
    assert np.allclose(matlab_round(np.array([0.5, -0.5, 1.4, 1.6])),
                       [1, -1, 1, 2])


def test_percentile():
    # MATLAB prctile of 1..10: P5=1, P50=5.5, P95=10
    x = np.arange(1, 11)
    assert percentile_matlab(x, 5) == 1.0
    assert percentile_matlab(x, 50) == 5.5
    assert percentile_matlab(x, 95) == 10.0


def test_fdr_bh():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205,
                  0.212, 0.216, 0.222, 0.251, 0.269, 0.275, 0.34])
    q = fdr_bh(p)
    m = len(p)
    s = np.sort(p)
    adj = np.minimum.accumulate((s * m / np.arange(1, m + 1))[::-1])[::-1]
    assert np.allclose(q, np.minimum(adj, 1))
    # monotone in p (weakly)
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12)
    # NaN handling
    assert np.isnan(fdr_bh(np.array([np.nan])))[0]


def test_grubbs_matches_isoutlier_semantics():
    rng = np.random.default_rng(3)
    x = rng.normal(size=200)
    x[[10, 150]] += 20          # two clear outliers (isoutlier finds both, iterated)
    mask = grubbs_outliers(x)
    assert set(np.flatnonzero(mask)) == {10, 150}
    # clean data: no outliers
    assert not grubbs_outliers(rng.normal(size=200)).any()
    # n < 3: no outliers
    assert not grubbs_outliers(np.array([1.0, 2.0])).any()


def test_valid_beats():
    # 1-based beats; window [-10 10] samples, half=5, n_pts=100
    b = np.array([1, 20, 90, 95])
    v = valid_beats(b, (-10, 10), 5, 100, None)
    # b=1: a=1-10-5=-14 <1 invalid; 20 ok; 90: z=90+10+5=105>100 invalid; 95 invalid
    assert list(np.flatnonzero(v)) == [1]
    # boundary
    bnd = [30]
    v = valid_beats(np.array([20, 22]), (-10, 10), 5, 100, bnd)
    # beat 22: window 12..37 contains 30 -> invalid; beat 20: 10..36 -> invalid too (30>10, 30<36)
    assert not v.any()
    # beat 25: window 10..35 contains 30 -> invalid
    assert not valid_beats(np.array([25]), (-10, 10), 5, 100, [30]).any()
    v = valid_beats(np.array([45]), (-10, 10), 5, 100, [30])
    assert v[0]          # window 30..50: bnd 30 not > 30 -> kept


def test_tf_stats_on_sinusoid():
    fs = 100.0
    f = 10.0
    t = np.arange(0, 2, 1 / fs)
    co = np.exp(2j * np.pi * f * t)          # unit sinusoid
    beats = np.arange(21, 100)               # inside
    t_idx = np.arange(1, 11)                 # arbitrary
    pw, pc = _tf_stats(co, beats, t_idx, co.size)
    assert np.allclose(pw, 1.0)
    # constant phase -> PPC = 1 (|sum z|^2 = N^2 -> 1 exactly)
    co1 = np.ones(t.size, dtype=complex)
    _, pc1 = _tf_stats(co1, beats, t_idx, co1.size)
    assert np.allclose(pc1, 1.0)
    # random phases -> PPC ~ 0 (unbiased)
    rng = np.random.default_rng(0)
    vals = []
    for _ in range(20):
        co2 = np.exp(1j * rng.uniform(0, 2 * np.pi, co.size))
        _, pc2 = _tf_stats(co2, beats, t_idx, co.size)
        vals.append(pc2.mean())
    assert abs(np.mean(vals)) < 0.05


def test_morlet_kernel_normalisation():
    fs = 250.0
    kernels, half = morlet_kernels_bb([4.0], fs, cycles=5.0)
    w = kernels[0]
    assert np.isclose(np.sum(np.abs(w)), 2.0)     # x2 amplitude normalisation
    # peak amplitude of |coef| for a sinusoid of amplitude A is A (test offline)
    # centre of the kernel is at the middle sample
    assert w.size % 2 == 1


def test_surrogate_rigid_min_shift():
    fs = 250.0
    beats = np.arange(100, 1000, 250)             # median IBI 250 samples (1 s)
    trains, shifts = surrogate_trains(beats, n_surr=20, mode='rigid', fs=fs,
                                      lims=(-75, 150), half=0, n_pts=4000)
    assert len(trains) == 20
    min_abs = int(matlab_round(0.25 * 250))       # 62 samples
    assert np.all(np.abs(shifts) >= min_abs)


def test_surrogate_shuffle_length():
    fs = 250.0
    all_beats = np.round(np.sort(np.cumsum(
        np.concatenate([[100], rng_int(60, 90, 300)])))).astype(int)
    beats = all_beats[valid_beats(all_beats, (-75, 150), 0, all_beats.max(), None)]
    trains, _ = surrogate_trains(beats, n_surr=5, mode='shuffle', fs=fs,
                                 all_beats=all_beats, lims=(-75, 150), half=0,
                                 n_pts=all_beats.max())
    for t in trains:
        assert len(t) <= len(beats)
        assert len(t) >= 10


def rng_int(lo, hi, n):
    return np.random.default_rng(1).integers(lo, hi, n)


def test_compute_hep_tf_end_to_end():
    """Synthetic: 10 Hz ring at 8 Hz bursts near R-peaks -> HRSP picks 8 Hz."""
    fs = 250.0
    n_pts = 60000
    t = np.arange(n_pts) / fs
    rng = np.random.default_rng(42)
    x = 0.5 * np.sin(2 * np.pi * 10 * t) + 0.05 * rng.normal(size=n_pts)
    # bursts at 8 Hz time-locked to beats
    beats = np.arange(2000, n_pts - 2000, 700, dtype=float)   # ~2.8 s apart
    for b in beats:
        i0 = int(b) - 100
        x[i0:i0 + 200] += 2.0 * np.sin(2 * np.pi * 8 * np.arange(200) / fs)
    tf, surr = compute_hep_tf(x, fs, beats, (-300, 600), freqs=np.arange(4, 21),
                              n_surr=0)
    # HEP at latency 0 has the burst average ~0 (sin over full cycles) but pow exists
    assert tf['hep'].shape == (1, len(tf['hep_times']))
    assert tf['hrsp'].shape == (1, 17, tf['times'].size)
    # surrogate comparison runs and produces finite p in [1/(n+1), 1]
    tf2, surr2 = compute_hep_tf(x, fs, beats, (-300, 600), freqs=np.arange(4, 13),
                                n_surr=20, seed=3)
    assert surr2['nSurr'] == 20
    assert np.nanmax(surr2['hrsp']['p_fdr'] <= 1.0)
    assert (surr2['hep']['p_fdr'] <= 1.0).all()
    # the 8 Hz burst sits at +0..+800 ms; HRSP vs shuffled nulls should show it
    # (weak assert: some point in the 8 Hz row significant)
    i_f = list(tf2['freqs']).index(8.0)
    assert (surr2['hrsp']['p_fdr'][0, i_f] < 0.05).any()