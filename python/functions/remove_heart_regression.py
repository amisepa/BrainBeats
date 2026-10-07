"""remove_heart_regression - Remove the cardiac field artifact from EEG by
regressing the (lagged) ECG out of each channel.

Python port of BrainBeats functions/remove_heart_regression.m (v1.6).
The ECG channel(s), filtered like the EEG and shifted by -20 to +20 ms
(one copy per sample), are fitted to each EEG channel by least squares over
the samples without large EEG artifacts (GFP within 6 robust SD of its
median); the fitted part is subtracted. No ICA needed; requires an ECG.

Returns (cleaned data, info) with info: lags_ms, samples_fitted (fraction in
%), cfa_before, cfa_after (heart-locked EEG amplitude, -50..100 ms, from
cfa_amplitude.m), variance_removed (%).

Copyright (C) Cedric Cannard, 2026 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import numpy as np

from .matlab_utils import matlab_round

__all__ = ['remove_heart_regression', 'cfa_amplitude']


def cfa_amplitude(D, beats, fs):
    """Heart-locked EEG amplitude (cfa_amplitude.m): RMS over channels of the
    heartbeat-locked average from -50 to 100 ms, relative to its mean from
    -200 to -100 ms. Beats are 1-based MATLAB sample indices."""
    D = np.asarray(D, float)
    w = np.arange(matlab_round(-0.2 * fs), matlab_round(0.4 * fs) + 1,
                  dtype=np.int64)
    beats = np.asarray(beats, np.int64)
    keep = (beats + w[0] > 0) & (beats + w[-1] <= D.shape[1])   # MATLAB 1-based
    beats = beats[keep]
    if beats.size < 10:
        return float('nan')
    pos = (beats[:, None] + w[None, :]) - 1                     # 0-based
    avg = D[:, pos].mean(axis=1)                                # nChan x nT
    t = w / fs
    avg = avg - avg[:, t < -0.1].mean(axis=1, keepdims=True)
    m = (t >= -0.05) & (t <= 0.1)
    return float(np.sqrt((avg[:, m] ** 2).mean()))


def remove_heart_regression(EEG, CARDIO, beats, params, *, filter_fn=None):
    """EEG/CARDIO: eegprep datasets (EEG channels only / ECG channels only).
    filter_fn(EEGorDict, lo, hi, causal) -> filtered dataset; defaults to
    eegprep.pop_eegfiltnew when available.

    Returns (EEG, info) -- EEG modified in place (its 'data' replaced).
    """
    import numbers
    max_lag = 20                                            # ms (remove_heart_regression.m:40)
    hp = params.get('highpass') or 0.5
    lp = params.get('lowpass') or 30.0
    causal = str(params.get('filttype', '')).lower() == 'causal'
    if filter_fn is None:
        import eegprep
        def filter_fn(eeg, lo, hi, causal):
            eeg = eegprep.pop_eegfiltnew(dict(eeg), locutoff=lo, minphase=causal)
            eeg = eeg[0] if isinstance(eeg, tuple) else eeg
            eeg = eegprep.pop_eegfiltnew(dict(eeg), hicutoff=hi, minphase=causal)
            return eeg[0] if isinstance(eeg, tuple) else eeg
    ECGd = filter_fn(dict(CARDIO), hp, lp, causal)
    ECGd = ECGd[0] if isinstance(ECGd, tuple) else ECGd
    n = min(int(EEG['pnts']), int(ECGd['pnts']))
    X = np.asarray(EEG['data'], float)[:, :n]
    ecg = np.asarray(ECGd['data'], float)[:, :n]
    fs = float(EEG['srate'])

    L = int(matlab_round(max_lag / 1000.0 * fs))
    n_reg = ecg.shape[0] * (2 * L + 1)
    R = np.empty((n_reg, n))
    k = 0
    for c in range(ecg.shape[0]):
        for lag in range(-L, L + 1):
            R[k] = np.roll(ecg[c], lag)
            k += 1
    R -= R.mean(axis=1, keepdims=True)

    gfp = X.std(axis=0, ddof=1)
    from scipy.stats import median_abs_deviation
    # MATLAB: mad(gfp,1) is the median absolute deviation; x1.4826 ~= robust SD
    ok = gfp < np.median(gfp) + 6.0 * median_abs_deviation(gfp, scale='normal')
    info = {'lags_ms': np.arange(-L, L + 1) / fs * 1000.0,
            'samples_fitted': float(100.0 * ok.mean())}
    Xok = X[:, ok]
    Rok = R[:, ok]
    B = (Xok @ Rok.T) @ np.linalg.pinv(Rok @ Rok.T)
    Y = X - B @ R
    # report on the artifact-free heartbeats (-200..400 ms clean)
    w = np.arange(matlab_round(-0.2 * fs), matlab_round(0.4 * fs) + 1,
                  dtype=np.int64)
    b = np.asarray(beats, np.int64)
    b = b[(b + w[0] >= 1) & (b + w[-1] <= n)]                  # 1-based
    b = b[np.asarray([bool(ok[(b_ + w) - 1].all()) for b_ in b])]
    info['cfa_before'] = cfa_amplitude(X, b, fs)
    info['cfa_after'] = cfa_amplitude(Y, b, fs)
    info['variance_removed'] = float(100.0 * (1 - (Y[:, ok] ** 2).sum() /
                                              (X[:, ok] ** 2).sum()))
    EEG['data'] = Y
    return EEG, info