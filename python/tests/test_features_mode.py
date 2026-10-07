"""Tests for the features mode (HRV + EEG features): MATLAB-parity of the
helpers and end-to-end sanity of the drivers against analytic expectations.

Run with the paa venv:  pytest test_features_mode.py
"""
import math
import sys

import numpy as np
import pytest

sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats-py/python')

from functions.compute_asymmetry import compute_asymmetry            # noqa: E402
from functions.compute_fe import compute_fe                          # noqa: E402
from functions.compute_psd import compute_psd, nextpow2              # noqa: E402
from functions.extract_features import extract_features              # noqa: E402
from functions.fractal_volatility import fractal_volatility          # noqa: E402
from functions.get_eeg_features import get_eeg_features              # noqa: E402
from functions.get_hrv_features import get_hrv_features              # noqa: E402
from functions.resample_NN import resample_NN                        # noqa: E402
from functions.params import default_params                          # noqa: E402
from functions.resting_iaf import resting_iaf, get_freqBounds        # noqa: E402
from functions.run_features import run_features                      # noqa: E402


# ---------------------------------------------------------------- nextpow2
def test_nextpow2():
    assert nextpow2(500) == 9 and 2 ** 9 >= 500 and 2 ** 8 < 500
    assert nextpow2(512) == 9


# ---------------------------------------------------------------- compute_psd
def test_compute_psd_peak():
    fs = 250
    x = np.sin(2 * np.pi * 10 * np.arange(2000) / fs)
    pwr, pwr_db, f = compute_psd(x, 500, 'hamming', 50, None, fs, [1, 40], 'psd')
    assert abs(f[np.argmax(pwr)] - 10) < 1.5
    assert np.allclose(pwr_db, 10 * np.log10(np.maximum(pwr, 2.2e-16)))
    # 'power' scaling keeps band integrals comparable to density x ENB
    p2, _, f2 = compute_psd(x, 500, 'hamming', 50, None, fs, [1, 40], 'power')
    assert p2.shape == pwr.shape


# ---------------------------------------------------------------- resample_NN
def test_resample_nn():
    t = np.arange(0, 20.0, 1.0)
    nn = np.linspace(1.0, 1.4, 20)
    y, tt = resample_NN(t, nn, 4, 'spline')
    assert tt[0] == 0 and len(tt) == 77   # 0:0.25:19 grid (MATLAB colon)
    # exact on the knots
    assert np.allclose(y[::4], nn, atol=1e-6)
    with pytest.raises(ValueError):
        resample_NN(t, nn, 4, 'nope')


# ---------------------------------------------------------------- compute_fe
def test_compute_fe_orders():
    rng = np.random.default_rng(42)
    t = np.arange(5000) / 250.
    fe_sin, p_sin = compute_fe(np.sin(2 * np.pi * 0.1 * t), 2, 0.15, 2, 1)
    fe_noise, p_noise = compute_fe(rng.normal(size=5000), 2, 0.15, 2, 1)
    assert 0 <= fe_sin < 0.1
    assert fe_noise > 1.0
    assert p_sin[0] > p_sin[1] and p_noise[0] > p_noise[1]   # p decreases with m
    fe_nan, p_nan = compute_fe(np.full(100, np.nan))
    assert math.isnan(fe_nan)


# ---------------------------------------------------------- fractal_volatility
def test_fractal_dimensions():
    rng = np.random.default_rng(1)
    fd_noise, sd1 = fractal_volatility(rng.normal(size=2000))
    fd_smooth, sd2 = fractal_volatility(np.sin(2 * np.pi * np.arange(4000) / 200.))
    assert 1.5 < fd_noise < 2.0
    assert 1.0 <= fd_smooth < 1.4
    assert sd1 > 0 and sd2 > 0


# --------------------------------------------------------------- compute_asymmetry
def test_asymmetry_basic():
    chanlocs = [{'labels': 'F3', 'X': 0.0, 'Y': 0.68, 'Z': 0.4},
                {'labels': 'F4', 'X': 0.0, 'Y': -0.68, 'Z': 0.4}]
    asy, labels, nums = compute_asymmetry(np.array([6.7, 0.6]), False, chanlocs)
    assert labels == ['F3 F4']
    assert nums.tolist() == [[0, 1]]
    assert abs(asy[0] - (math.log(6.7) - math.log(0.6))) < 1e-6
    # normalized variant: divide by each channel's total power first
    asy2, _, _ = compute_asymmetry(np.array([6.7, 0.6]), True, chanlocs,
                                   False, np.array([10.0, 2.0]))
    assert abs(asy2[0] - (math.log(0.67) - math.log(0.30))) < 1e-6


# ------------------------------------------------------------ get_hrv_features
def _synthetic_nn(seed=11, n=1800):
    rng = np.random.default_rng(seed)
    nn = 1.0 + 0.02 * rng.standard_normal(n)
    return nn, np.cumsum(nn) - nn[0]


def test_hrv_time_domain_matches_hand():
    nn, nn_t = _synthetic_nn()
    hrv, p = get_hrv_features(nn, nn_t, {'hrv_time': True})
    assert hrv['time']['SDNN'] == round(float(np.std(nn * 1000, ddof=1)), 1)
    assert hrv['time']['RMSSD'] == round(
        float(np.sqrt(np.mean(np.diff(nn * 1000) ** 2))), 1)


def test_hrv_frequency_and_nonlinear():
    nn, nn_t = _synthetic_nn()
    hrv, p = get_hrv_features(nn, nn_t, {'hrv_time': True, 'hrv_frequency': True,
                                         'hrv_nonlinear': True})
    fr = hrv['frequency']
    assert set(fr) >= {'vlf', 'lf', 'hf', 'lfhf'}      # ULF skipped (30-min file)
    assert 0 < fr['lf'] < fr['hf']                     # 0.02-SD NN: HF dominates
    nl = hrv['nonlinear']
    assert abs(nl['Poincare']['SD1'] - 20.0) < 1.5     # 20 ms RMSSD-like noise
    assert 0.5 < nl['FD'] < 2.0 and 0.5 < nl['FE'] < 2.0
    assert np.isfinite(nl['PRSA_AC']) and np.isfinite(nl['PRSA_DC'])


def test_hrv_poincare_analytic():
    """Poincare: SD1^2 = var(diff)/2; SD2^2 = 2*var - var(diff)/2 (1 s IBI)."""
    nn, nn_t = _synthetic_nn(seed=3)
    hrv, _ = get_hrv_features(nn, nn_t, {'hrv_nonlinear': True})
    d = np.diff(nn)
    sd1 = math.sqrt(np.var(d, ddof=1) / 2) * 1000
    assert abs(hrv['nonlinear']['Poincare']['SD1'] - round(sd1, 3)) < 0.01


# ---------------------------------------------------------- get_eeg_features
def test_eeg_features_core():
    rng = np.random.default_rng(5)
    fs = 250
    n = fs * 300
    t = np.arange(n) / fs
    X = 5 * rng.standard_normal((3, n))
    X[0] += 8 * np.sin(2 * np.pi * 10 * t)          # alpha on F3
    chanlocs = [{'labels': 'F3', 'X': 0.0, 'Y': 0.68, 'Z': 0.4},
                {'labels': 'F4', 'X': 0.0, 'Y': -0.68, 'Z': 0.4},
                {'labels': 'Cz', 'X': 0.0, 'Y': 0.02, 'Z': 0.95}]
    p = {'fs': fs, 'chanlocs': chanlocs, 'eeg_time': True,
         'eeg_frequency': True, 'eeg_nonlinear': False,
         'eeg_frange': [1, 40], 'eeg_norm': 0}
    feats, p2 = get_eeg_features(X, p)
    tm = feats['time']
    assert np.allclose(tm['rms'], np.sqrt((X ** 2).mean(axis=1)), rtol=1e-9)
    fr = feats['frequency']
    # 10-Hz peak channel (F3) has by far the largest alpha power
    assert fr['alpha'][0] > 5 * fr['alpha'][1]
    # qEEG: rel powers sum to ~1 per channel
    # band gaps (3-4, 7-8, 13-13 shared) -> rel sums slightly below 1
    tot = sum(feats['qeeg'][k] for k in ('delta_rel', 'theta_rel', 'alpha_rel',
                                         'beta_rel', 'gamma_rel'))
    assert np.all((tot > 0.90) & (tot <= 1.0 + 1e-9))
    # peak-alpha of the injected channel = near 10 Hz (MATLAB argmax in [7 13])
    assert abs(feats['qeeg']['iaf'][0] - 10) < 1.5
    # alpha asymmetry: F3 > F4 -> positive asy
    idx = {lab: i for i, lab in enumerate(
        feats['frequency']['asymmetry_pairs_labels'])}
    iAsy = list(feats['frequency']['asymmetry_pairs_labels']).index('F3 F4')
    assert feats['frequency']['asymmetry'][iAsy] > 1.0
    assert np.isfinite(feats['qeeg']['median']).all()
    assert np.isfinite(feats['qeeg']['sef90']).all()


def test_eeg_individualized_alpha_bounds():
    fs = 100
    n = fs * 180
    t = np.arange(n) / fs
    X = 5 * np.random.default_rng(7).standard_normal((1, n)) \
        + 6 * np.sin(2 * np.pi * 10 * t)
    p = {'fs': fs, 'chanlocs': [{'labels': 'Cz', 'X': 0, 'Y': 0.02, 'Z': 0.95}],
         'eeg_time': False, 'eeg_frequency': True, 'eeg_nonlinear': False,
         'eeg_frange': [1, 40], 'eeg_norm': 0,
         'eeg_freqbounds': 'individualized'}
    feats, _ = get_eeg_features(X, p)
    bands = feats['frequency']['bands']
    lo, hi = bands[2]
    assert 7 <= lo <= 10 and 10 <= hi <= 13            # around the 10-Hz peak


# ------------------------------------------------------------ restingIAF core
def test_resting_iaf_recovers_alpha():
    fs = 250
    n = fs * 300
    t = np.arange(n) / fs
    rng = np.random.default_rng(5)
    X = 8 * np.sin(2 * np.pi * 10 * t) + 5 * rng.standard_normal(n)
    pSum, pChans = resting_iaf(X[None, :], [1, 30], [7, 14], 11, 5, fs=fs, cmin=1)
    assert abs(pSum['cog'] - 10.0) < 0.6
    # get_freqBounds on the same spectrum: bounds bracket 10 Hz
    from functions.compute_psd import compute_psd
    pwr, _, f = compute_psd(X, fs * 4, 'hamming', 50, None, fs, [1, 30], 'psd')
    bounds, peak = get_freqBounds(pwr, f, fs, [7, 14], fs * 4, 1)
    assert np.isfinite(peak) and abs(peak - 10) < 1.0


# ------------------------------------------------------------- extract_features
def test_extract_features_tables():
    Features = {
        'HRV': {'time': {'SDNN': 19.9, 'RMSSD': 28.3, 'pNN50': 7.2},
                'frequency': {'ulf': 1.0, 'vlf': 2.0, 'lf': 3.0, 'hf': 4.0,
                              'lfhf': 0.75},
                'nonlinear': {'Poincare': {'SD1': 20.0, 'SD2': 19.9,
                                           'SD1SD2': 1.006},
                              'PRSA_AC': -5.7, 'PRSA_DC': 5.9, 'FE': 1.28,
                              'FD': 1.8}},
        'EEG': {'time': {'rms': np.array([5.0, 7.0])},
                'frequency': {'delta': np.array([1., 2.]),
                              'IAF': np.array([10.1, 10.2]), 'IAF_mean': 10.15,
                              'asymmetry': np.array([2.41]),
                              'asymmetry_pairs_labels': ['F3 F4']},
                'qeeg': {'total_abs': np.array([7.6, 7.8])},
                'nonlinear': {'FD': np.array([1.8, 1.9])}},
    }
    hrv, eeg = extract_features(Features)
    assert hrv['SDNN-HRV'][0] == 19.9 and 'HF-HRV' in hrv
    assert abs(hrv['LF/HF'][0] - 0.75) < 1e-9
    assert hrv['Poincaré: SD1'][0] == 20.0 and 'HRV-MFE_peak' not in hrv
    assert eeg['EEG-rms_1'][0] == 5.0 and eeg['EEG-rms_2'][0] == 7.0
    assert eeg['EEG-IAF_mean'][0] == 10.15
    assert abs(eeg['Asy (F3 F4)'][0] - 2.41) < 1e-9
    assert eeg['qEEG-total_abs_2'][0] == 7.8


# ------------------------------------------------------------- run_features
def test_run_features_endtoend():
    rng = np.random.default_rng(3)
    fs = 100
    n = fs * 180
    t = np.arange(n) / fs
    ecg = 50 * rng.standard_normal(n)
    eeg = 5 * rng.standard_normal((1, n)) + 4 * np.sin(2 * np.pi * 10 * t)
    # NN input directly (skip ECG peak detection dependency)
    nn = 1.0 + 0.02 * rng.standard_normal(120)
    nn_t = np.cumsum(nn) - nn[0]
    p = default_params(analysis='features', fs=fs,
                       chanlocs=[{'labels': 'Cz', 'X': 0, 'Y': 0.02, 'Z': .95}])
    p['hrv_time'] = p['hrv_frequency'] = p['hrv_nonlinear'] = True
    p['eeg_time'] = p['eeg_frequency'] = True
    p['eeg_nonlinear'] = True
    feats, p2 = run_features(ECG=ecg, NN=nn, NN_times=nn_t, EEG=eeg, params=p)
    assert 'HRV' in feats and 'EEG' in feats
    assert 'SDNN' in feats['HRV']['time']
    assert np.isfinite(feats['EEG']['frequency']['IAF']).all()
    assert np.isfinite(feats['EEG']['nonlinear']['FE']).all()
    # features tables round-trip
    hrv, eegTab = extract_features({'HRV': feats['HRV'], 'EEG': feats['EEG']})
    assert 'SDNN-HRV' in hrv and 'EEG-IAF' in eegTab