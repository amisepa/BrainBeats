"""Unit tests for the v1.6 modules: estimate_pat, estimate_pat_eeg,
csd_transform, cfa_amplitude."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))          # the python/ package root

from functions.estimate_pat import estimate_pat               # noqa: E402
from functions.estimate_pat_eeg import estimate_pat_eeg       # noqa: E402
from functions.csd_transform import (read_ced_chanlocs, get_gh,  # noqa: E402
                                     csd_matrix, csd_transform_data)
from functions.remove_heart_regression import cfa_amplitude   # noqa: E402

FS = 250.0


def _heart_train(n_sec=30.0, ibi_s=1.0, fs=FS, offset=0):
    beats = np.arange(offset, n_sec * fs, ibi_s * fs)
    return beats.astype(np.int64)


def test_estimate_pat_pairs_preceding_rpeak():
    rpk = _heart_train(n_sec=120, offset=100)      # R-peaks every 1 s in 120 s
    ppg = rpk + int(0.24 * FS)                     # each pulse 240 ms later
    med, pat, info = estimate_pat(rpk, ppg, FS)
    assert abs(med - 240.0) < 1.0
    assert info['nPaired'] == rpk.size             # all paired (240 in [50 600])
    assert np.allclose(pat[np.isfinite(pat)], 240.0, atol=1e-6)


def test_estimate_pat_range_filter():
    rpk = _heart_train(n_sec=120, offset=100)
    ppg = rpk + int(0.8 * FS)                      # 800 ms: outside [50 600]
    med, pat, info = estimate_pat(rpk, ppg, FS)
    assert np.isnan(med) == (info['nPaired'] == 0)
    assert info['nPaired'] == 0


def test_estimate_pat_eeg_recovers_cfa_delay():
    fs = FS
    beats = _heart_train(39.0, 1.0, fs, offset=int(1.5 * fs))
    n = int(40.0 * fs)
    t = np.arange(n) / fs
    x = np.zeros((3, n))
    delay = 0.35                                    # s (PAT to recover)
    # QRS-like kernel: sharp positive then slower negative (the cardiac field)
    tt = np.arange(-0.05, 0.15, 1 / fs)
    kernel = np.exp(-((tt - 0.0) ** 2) / (2 * 0.012 ** 2)) \
        - 0.3 * np.exp(-((tt - 0.09) ** 2) / (2 * 0.05 ** 2))
    for b in beats:
        i0 = b - int(delay * fs)            # QRS at the PULSE minus PAT
        i1 = min(i0 + kernel.size, x.shape[1])
        for ch in range(3):
            x[ch, i0:i1] += kernel[:i1 - i0] * (1.0 + 0.1 * ch)
    pat, info = estimate_pat_eeg(x, fs, beats)
    # the GFP peak should sit at the kernel onset shifted by -delay
    assert info['ratio'] >= 2
    assert pat == pytest.approx(delay * 1000, abs=60)   # within 60 ms


def test_get_gh_and_csd_matrix_shapes_and_symmetry():
    labels, theta, phi = read_ced_chanlocs(os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'functions', 'standard_1005.ced'))
    assert len(labels) > 100
    sub = slice(0, 19)
    G, H = get_gh(theta[sub], phi[sub], 4)
    assert G.shape == H.shape == (19, 19)
    assert np.allclose(G, G.T, atol=1e-10)
    assert np.allclose(H, H.T, atol=1e-10)
    C = csd_matrix(19, G, H)
    assert C.shape == (19, 19)


def test_csd_transform_data_applies_linear_transform():
    labels, theta, phi = read_ced_chanlocs(os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'functions', 'standard_1005.ced'))
    keep = ['Fz', 'Cz', 'Pz']
    idx = [labels.index(l) for l in keep]
    rng = np.random.default_rng(3)
    data = rng.normal(size=(3, 50))
    out, C = csd_transform_data(data, keep)
    assert out.shape == data.shape
    assert np.allclose(out, C @ data)
    # a spatially uniform field (all channels equal) maps to ~0 (reference-free)
    flat = np.ones((3, 50))
    out_f, _ = csd_transform_data(flat, keep)
    assert np.abs(out_f).max() < 1e-6


def test_cfa_amplitude_rms_formula():
    fs = FS
    n = int(14 * fs)
    beats = _heart_train(13.5, 1.0, fs, offset=int(1.0 * fs))
    assert beats.size >= 10                  # cfa_amplitude needs 10 heartbeats
    x = np.zeros((2, n))
    tt = np.arange(-0.2, 0.4, 1 / fs)
    kernel = np.exp(-((tt) ** 2) / (2 * 0.02 ** 2))     # peak at t=0, 1 uV tall
    for b in beats:
        i0 = b - int(0.2 * fs)
        i1 = min(i0 + tt.size, x.shape[1])
        for ch in range(2):
            x[ch, i0:i1] += kernel[:i1 - i0] * (1.0 + ch)
    a = cfa_amplitude(x, beats, fs)
    # replay the exact cfa_amplitude formula on the construction: grand
    # average over beats (all windows fit), baseline -200..-100 ms (~0, the
    # kernel is Gaussian sigma=20 ms), then RMS over channels and the
    # -50..100 ms latencies
    off = np.arange(int(round(-0.2 * fs)), int(round(0.4 * fs)) + 1)
    tw = off / fs
    pos = (beats[:, None] + off[None, :]) - 1
    ga = x[:, pos].mean(axis=1)                       # channels x lags
    base = ga[:, tw < -0.1].mean(axis=1, keepdims=True)
    mm = (tw >= -0.05) & (tw <= 0.1)
    expected = np.sqrt(np.mean((ga - base)[:, mm] ** 2))
    assert a == pytest.approx(expected, rel=1e-9)


def test_cfa_amplitude_needs_10_beats():
    fs = FS
    x = np.zeros((2, int(3 * fs)))
    assert np.isnan(cfa_amplitude(x, np.arange(1, 5), fs))


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-q', '-p', 'no:cacheprovider']))