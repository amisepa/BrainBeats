"""get_hrv_features - Heart-rate variability features from an NN series.

Python port of BrainBeats functions/get_hrv_features.m. Time domain
(SDNN, RMSSD, pNN50), frequency domain (ULF/VLF/LF/HF band powers in
sliding windows of 5 cycles of the band's lowest frequency, Task Force
1996; Lomb-Scargle via the validated astropy implementation, Welch/FFT on
7-Hz resampled NN) and nonlinear domain (Poincare SD1/SD2, fuzzy entropy,
fractal dimension, PRSA acceleration/deceleration capacity).

Copyright (C) Cedric Cannard, 2023 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import numpy as np
from astropy.timeseries import LombScargle

from .compute_fe import compute_fe
from .compute_psd import nextpow2
from .compute_psd import compute_psd as _compute_psd
from .fractal_volatility import fractal_volatility
from .resample_NN import resample_NN

__all__ = ['get_hrv_features']


def _round_sig(x, nd):
    """MATLAB round(x, n, 'significant')."""
    if not np.isfinite(x):
        return float(x)
    if x == 0:
        return 0.0
    from math import floor, log10 as lg
    return round(x, -int(floor(lg(abs(x)))) + (nd - 1))


def get_hrv_features(NN, NN_times, params):
    """NN (s), NN_times (s). Reads params['hrv_time'/'hrv_frequency'/
    'hrv_nonlinear'], optional hrv_spec ('LombScargle_norm' default),
    hrv_norm (false), hrv_overlap (0.25). Returns (HRV dict, params)."""
    NN = np.asarray(NN, float).ravel()
    NN_times = np.asarray(NN_times, float).ravel()
    HRV = {}

    # ---- Time domain ----------------------------------------------------
    if params.get('hrv_time'):
        print('Extracting HRV features in the time domain...')
        nn_ms = NN * 1000.0
        HRV['time'] = {
            'SDNN': round(float(np.std(nn_ms, ddof=1)), 1),
            'RMSSD': round(float(np.sqrt(np.mean(np.diff(nn_ms) ** 2))), 1),
            'pNN50': round(float(100.0 * np.mean(np.abs(np.diff(NN)) >= 0.050)), 1),
        }

    # ---- Frequency domain -----------------------------------------------
    if params.get('hrv_frequency'):
        print('Extracting HRV features in the frequency domain...')
        hrv_spec = str(params.get('hrv_spec', 'LombScargle_norm'))
        if hrv_spec not in ('LombScargle_norm', 'LombScargle', 'welch', 'fft'):
            raise ValueError(f"get_hrv_features: unknown hrv_spec '{hrv_spec}'.")
        norm = bool(params.get('hrv_norm', False))
        overlap = float(params.get('hrv_overlap', 0.25))
        params['hrv_spec'] = hrv_spec
        params['hrv_norm'] = norm
        params['hrv_overlap'] = overlap

        pwrScale = 1.0 if hrv_spec == 'LombScargle_norm' else 1e6
        bands = np.array([[0.0, 0.003], [0.003, 0.04], [0.04, 0.15], [0.15, 0.40]])
        bandNames = ['ULF', 'VLF', 'LF', 'HF']
        params['hrv_band_freqs'] = bands
        params['hrv_band_names'] = bandNames
        minLength = np.ceil([86400.0, 5 / 0.003, 5 / 0.04, 5 / 0.15])

        HRV['frequency'] = {}
        PWR, PWR_freqs = {}, {}
        anyBand = False
        for iBand in range(4):
            lo, hi = bands[iBand]
            if NN_times[-1] < minLength[iBand]:
                print(f'WARNING: file too short for reliable {bandNames[iBand]} '
                      f'power ({minLength[iBand] / 60:g} min needed): skipped.')
                continue
            anyBand = True
            winLength = float(minLength[iBand])
            stepSize = float(np.floor(winLength * (1 - overlap)))
            nWindows = int(np.floor((NN_times[-1] - winLength) / stepSize)) + 1
            per_win = []
            for iWin in range(nWindows):
                start_idx = iWin * stepSize + 1.0
                end_idx = start_idx + winLength - 1.0
                win_idx = (NN_times >= start_idx) & (NN_times <= end_idx)
                if NN[win_idx].sum() < 0.85 * winLength:
                    print('WARNING: this HRV window contains a gap > 15% of the '
                          'band window (likely from RR artifact cleaning).')
                n_in = int(win_idx.sum())
                if n_in < 2:
                    continue
                if hrv_spec in ('LombScargle_norm', 'LombScargle'):
                    nfft = 2 ** nextpow2(n_in)
                    freqs = np.arange(lo, hi + 1.0 / nfft, 1.0 / nfft)
                    ls = LombScargle(NN_times[win_idx], NN[win_idx],
                                     normalization='standard',
                                     center_data=bool(hrv_spec == 'LombScargle_norm'))
                    pwr = ls.power(freqs)          # unitless, 0..1 normalization
                    if hrv_spec == 'LombScargle_norm':
                        # MATLAB plomb 'normalized': periodogram divided by the
                        # series variance (variance-normalized PSD)
                        pwr = pwr / 2.0            # astropy: sum(pwr)*df ~ var/2
                else:
                    resamp_freq = 7
                    NN_resamp, _ = resample_NN(NN_times[win_idx], NN[win_idx],
                                               resamp_freq, 'cub')
                    NN_resamp = NN_resamp - NN_resamp.mean()
                    if hrv_spec == 'welch':
                        welchWin = int(min(NN_resamp.size,
                                           round(minLength[iBand] * resamp_freq)))
                        pwr, freqs, _ = _compute_psd(NN_resamp, welchWin,
                                                     'hamming', 50, None,
                                                     resamp_freq,
                                                     [bands[iBand][0],
                                                      bands[iBand][1]], 'psd')
                    else:                          # fft
                        nSamp = NN_resamp.size
                        nHalf = nSamp // 2 + 1
                        pwr_all = np.abs(np.fft.fft(NN_resamp)) ** 2 / (resamp_freq * nSamp)
                        pwr_all = pwr_all[:nHalf]
                        pwr_all[1:-1 + (nSamp % 2)] *= 2
                        freqs = np.arange(nHalf) * resamp_freq / nSamp
                keep = (freqs >= lo) & (freqs <= hi)
                freq_res = freqs[1] - freqs[0] if freqs.size > 1 else 1.0
                per_win.append(float(np.sum(pwr[keep]) * freq_res * pwrScale))
                PWR[bandNames[iBand]] = np.asarray(pwr[keep])
                PWR_freqs[bandNames[iBand]] = freqs[keep]
            if per_win:
                HRV['frequency'][bandNames[iBand].lower()] = np.asarray(per_win)

        if not anyBand:
            print('WARNING: recording too short for any HRV frequency band '
                  '(at least 34 s needed). Band powers set to NaN.')
            HRV['frequency'] = {'ulf': np.nan, 'vlf': np.nan, 'lf': np.nan,
                                'hf': np.nan}

        freq = HRV['frequency']
        if 'lf' in freq and 'hf' in freq:
            freq['lfhf'] = round(float(np.mean(freq['lf']) / np.mean(freq['hf'])), 2)
        if norm:
            try:
                freq['ttlpwr'] = float(np.sum([np.mean(freq[k]) for k in
                                               ('ulf', 'vlf', 'lf', 'hf')
                                               if k in freq]))
            except Exception:
                freq['ttlpwr'] = float(np.sum([np.mean(freq[k]) for k in
                                               ('vlf', 'lf', 'hf') if k in freq]))
            for k in ('ulf', 'vlf', 'lf', 'hf'):
                if k in freq:
                    freq[k] = float(np.mean(freq[k]) / freq['ttlpwr'])
        if PWR:
            freq['pwr'] = np.concatenate([PWR[k] for k in PWR])
            freq['pwr_freqs'] = np.concatenate([PWR_freqs[k] for k in PWR_freqs])
            freq['bands'] = bands
        rnd = (lambda x: _round_sig(x, 4)
               ) if (hrv_spec == 'LombScargle_norm' or norm) else (lambda x: round(x, 2))
        for k in ('ulf', 'vlf', 'lf', 'hf'):
            if k in freq and np.ndim(freq[k]) > 0:
                freq[k] = rnd(float(np.nanmean(freq[k])))
        HRV['frequency'] = freq

    # ---- Nonlinear domain ------------------------------------------------
    if params.get('hrv_nonlinear'):
        print('Extracting HRV features in the nonlinear domain (Poincare, '
              'fuzzy entropy, fractal dimension, PRSA)...')
        SDSD = float(np.std(np.diff(NN), ddof=1))
        SDRR = float(np.std(NN, ddof=1))
        SD1 = SDSD / np.sqrt(2.0)
        SD2 = np.sqrt(2.0 * SDRR ** 2 - 0.5 * SDSD ** 2)
        HRV['nonlinear'] = {
            'Poincare': {'SD1': round(SD1 * 1000.0, 3),
                         'SD2': round(SD2 * 1000.0, 3),
                         'SD1SD2': round(SD1 / SD2, 3)},
        }
        m, r, n, tau = 2, 0.15, 2, 1
        params['entropy_m'], params['entropy_r'] = m, r
        params['entropy_tau'], params['entropy_n'] = tau, n
        fe, _ = compute_fe(NN, m, r, n, tau)
        HRV['nonlinear']['FE'] = fe
        fd, _ = fractal_volatility(NN)
        HRV['nonlinear']['FD'] = fd
        print('Computing phase rectified signal averaging (PRSA)...')
        thresh = 20
        params['prsa_thresh'] = thresh
        lowAnchor = 1 - thresh / 100 - 0.0001
        highAnchor = 1 + thresh / 100
        drr_per = NN[1:] / NN[:-1]
        ac_anchor = np.flatnonzero((drr_per > lowAnchor) & (drr_per <= 0.9999)) + 1
        dc_anchor = np.flatnonzero((drr_per > 1) & (drr_per <= highAnchor)) + 1
        HRV['nonlinear']['PRSA_AC'] = _prsa_capacity(NN, ac_anchor)
        HRV['nonlinear']['PRSA_DC'] = _prsa_capacity(NN, dc_anchor)

    return HRV, params


def _prsa_capacity(NN, anchors):
    """Acceleration/deceleration capacity (ms); drop anchors without two
    beats before and one after: X = [X(-2) X(-1) X(0) X(1)] around each
    anchor, cap = (X(0) + X(1) - X(-1) - X(-2)) / 4."""
    anchors = np.asarray(anchors, int)
    # anchors point at 0-based index i; need i-2 >= 0 and i+1 <= n-1
    a_ok = anchors[(anchors - 2 >= 0) & (anchors + 1 <= NN.size - 1)]
    if a_ok.size == 0:
        return np.nan
    win = np.array([NN[a - 2:a + 2] for a in a_ok])       # (nAnchors, 4)
    X = win.mean(axis=0)                                  # [X(-2)..X(1)]
    return round(1000.0 * (X[2] + X[3] - X[1] - X[0]) / 4.0, 2)