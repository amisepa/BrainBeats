"""get_eeg_features - EEG features in the time, frequency and nonlinear
domains for each channel.

Python port of BrainBeats functions/get_eeg_features.m: time (rms, mode,
var, skewness, kurtosis, iqr); frequency (Welch PSD, conventional or
individualized band powers, qEEG absolute/relative powers and ratios, IAF
via the restingIAF method, alpha asymmetry); nonlinear (fractal dimension,
fuzzy entropy at m=2, r=.15, n=2, tau=1 on data resampled to ~90 Hz when
fs > 100 Hz).

Copyright (C) Cedric Cannard, 2023 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import iqr, kurtosis, skew

from .compute_asymmetry import compute_asymmetry
from .compute_fe import compute_fe
from .compute_psd import compute_psd
from .fractal_volatility import fractal_volatility
from .resting_iaf import get_freqBounds, resting_iaf
from .resample_NN import resample_NN  # noqa: F401  (parity helper)

__all__ = ['get_eeg_features']


def _mode_row(x):
    """MATLAB mode(signals, 2): most frequent value per row."""
    vals, counts = np.unique(x[~np.isnan(x)], return_counts=True)
    return vals[np.argmax(counts)] if vals.size else np.nan


def _resample_90(signals, fs):
    """MATLAB resample(sig, p, q) to 90 Hz via rat(new_fs/fs, 1e-4)."""
    from fractions import Fraction
    fr = Fraction(90.0 / float(fs)).limit_denominator(10 ** 9)
    # MATLAB rat with 1e-4 tolerance: find smallest p/q with |p/q - 90/fs| < 1e-4
    p_found, q_found = None, None
    for q in range(1, 5001):
        p = round(90.0 * q / fs)
        if p == 0:
            continue
        err = abs(p / q - 90.0 / fs)
        if err < 1e-4:
            p_found, q_found = p, q
            break
    if p_found is None:
        p_found, q_found = fr.numerator, fr.denominator
    from scipy.signal import resample as sp_resample
    n_out = int(np.ceil(signals.shape[1] * p_found / q_found))
    out = np.empty((signals.shape[0], n_out))
    for i in range(signals.shape[0]):
        out[i] = sp_resample(signals[i], n_out, window=('kaiser', 5.0))
    return out, fs * p_found / q_found


def get_eeg_features(signals, params):
    """signals: (channels x samples, uV). Reads params eeg_time/eeg_frequency/
    eeg_nonlinear (+fs, chanlocs; defaults eeg_frange [1 40], eeg_wintype
    'hamming', eeg_winlen 2, eeg_winoverlap 50, eeg_freqbounds
    'conventional', eeg_norm 1, asy_norm False). Returns (features dict,
    params)."""
    signals = np.asarray(signals, float)
    n_chan = signals.shape[0]
    fs = params.get('fs')
    if not fs:
        raise ValueError("get_eeg_features: params['fs'] is required")
    feats = {}

    # ---- Time domain ------------------------------------------------------
    if params.get('eeg_time'):
        print('Calculating time-domain EEG features...')
        feats['time'] = {
            'rms': np.sqrt((signals ** 2).mean(axis=1)),
            'mode': np.array([_mode_row(s) for s in signals]),
            'var': signals.var(axis=1, ddof=1),
            'skewness': skew(signals, axis=1, bias=False, nan_policy='omit'),
            'kurtosis': kurtosis(signals, axis=1, fisher=False, bias=False,
                                 nan_policy='omit'),
            'iqr': iqr(signals, axis=1),
        }

    # ---- Frequency domain -------------------------------------------------
    if params.get('eeg_frequency'):
        fRange = params.get('eeg_frange') or [1, 40]
        params.setdefault('eeg_frange', fRange)
        wintype = params.get('eeg_wintype') or 'hamming'
        params.setdefault('eeg_wintype', wintype)
        winlen = params.get('eeg_winlen') or 2
        params.setdefault('eeg_winlen', winlen)
        overlap = params.get('eeg_winoverlap') or 50
        params.setdefault('eeg_winoverlap', overlap)
        freqbounds = str(params.get('eeg_freqbounds', 'conventional')).lower()
        params.setdefault('eeg_freqbounds', freqbounds)
        eeg_norm = params.get('eeg_norm', 1)
        params.setdefault('eeg_norm', eeg_norm)
        asy_norm = params.get('asy_norm', False)
        params.setdefault('asy_norm', asy_norm)

        _, _, f = compute_psd(signals[0], fs * winlen, wintype, overlap,
                              None, fs, fRange, 'psd')
        f = np.asarray(f)
        bands = np.array([[f[0], 3.0], [4.0, 7.0], [8.0, 13.0],
                          [13.0, 30.0], [30.0, fRange[1]]])

        if freqbounds == 'individualized':
            print('Estimating individualized alpha band bounds...')
            pwr_all, _, _ = compute_psd(signals, fs * winlen, wintype, overlap,
                                        None, fs, fRange, 'psd')
            bounds = []
            for iC in range(n_chan):
                try:
                    b, _ = get_freqBounds(pwr_all[iC], f, fs, [7, 14],
                                          fs * winlen, 1)
                    if np.all(np.isfinite(b)):
                        bounds.append(b)
                except Exception:
                    pass
            if bounds:
                alphaBounds = np.median(np.asarray(bounds), axis=0)
                if (np.all(np.isfinite(alphaBounds))
                        and alphaBounds[0] > bands[1, 0]
                        and alphaBounds[1] < bands[3, 1]):
                    bands[1, 1] = alphaBounds[0]
                    bands[2, :] = alphaBounds
                    bands[3, 0] = alphaBounds[1]
                    print(f'Individualized alpha band: {alphaBounds[0]:.2f}-'
                          f'{alphaBounds[1]:.2f} Hz')
                else:
                    print('WARNING: no alpha peak detected to individualize the '
                          "frequency bands. Using conventional bands.")

        df = float(np.mean(np.diff(f)))
        PWR = np.empty((n_chan, f.size))
        band_pwr = {b: np.empty(n_chan) for b in
                    ('delta', 'theta', 'alpha', 'beta', 'gamma')}
        qeeg = {k: np.empty(n_chan) for k in
                ('delta_abs', 'theta_abs', 'alpha_abs', 'beta_abs', 'gamma_abs',
                 'total_abs', 'delta_rel', 'theta_rel', 'alpha_rel', 'beta_rel',
                 'gamma_rel', 'alpha_theta', 'theta_beta', 'alpha_beta',
                 'alpha_tplusb', 'iaf', 'median', 'sef90')}
        n_idx = [f >= bands[0, 0], f <= bands[1, 1], np.zeros(f.size, bool)]
        idxD = (f >= bands[0, 0]) & (f <= bands[0, 1])
        idxT = (f >= bands[1, 0]) & (f <= bands[1, 1])
        idxA = (f >= bands[2, 0]) & (f <= bands[2, 1])
        idxB = (f >= bands[3, 0]) & (f <= bands[3, 1])
        idxG = (f >= bands[4, 0]) & (f <= bands[4, 1])
        print('Calculating band-power on each EEG channel:')
        for iC in range(n_chan):
            pwr, pwr_db, _ = compute_psd(signals[iC], fs * winlen, wintype,
                                         overlap, None, fs, fRange, 'psd')
            PWR[iC] = pwr
            if eeg_norm == 0:
                bp = [pwr[m].mean() for m in (idxD, idxT, idxA, idxB, idxG)]
            elif eeg_norm == 1:
                bp = [pwr_db[m].mean() for m in (idxD, idxT, idxA, idxB, idxG)]
            else:
                bp = [pwr[m].mean() for m in (idxD, idxT, idxA, idxB, idxG)]
                bp = np.asarray(bp) / pwr.sum()
            for key, v in zip(('delta', 'theta', 'alpha', 'beta', 'gamma'), bp):
                band_pwr[key][iC] = v
            tot = float(pwr.sum() * df)
            qeeg['total_abs'][iC] = tot
            absP = [float(pwr[m].sum() * df) for m in (idxD, idxT, idxA, idxB, idxG)]
            for key, v in zip(('delta_abs', 'theta_abs', 'alpha_abs', 'beta_abs',
                               'gamma_abs'), absP):
                qeeg[key][iC] = v
            den = max(tot, np.finfo(float).eps)
            for key, v in zip(('delta_rel', 'theta_rel', 'alpha_rel', 'beta_rel',
                               'gamma_rel'), absP):
                qeeg[key][iC] = v / den
            qeeg['alpha_theta'][iC] = absP[2] / max(absP[1], np.finfo(float).eps)
            qeeg['theta_beta'][iC] = absP[1] / max(absP[3], np.finfo(float).eps)
            qeeg['alpha_beta'][iC] = absP[2] / max(absP[3], np.finfo(float).eps)
            qeeg['alpha_tplusb'][iC] = (absP[2] / max(absP[1] + absP[3],
                                                      np.finfo(float).eps))
            aSearch = (f >= 7) & (f <= 13)
            qeeg['iaf'][iC] = (f[aSearch][np.argmax(pwr[aSearch])]
                               if aSearch.any() else np.nan)
            cs = np.cumsum(pwr) * df
            if cs[-1] > 0:
                qeeg['median'][iC] = np.interp(0.5 * cs[-1], cs, f)
                qeeg['sef90'][iC] = np.interp(0.9 * cs[-1], cs, f)
            else:
                qeeg['median'][iC] = qeeg['sef90'][iC] = np.nan
        feats['frequency'] = {'freqs': f, 'pwr': PWR, 'bands': bands}
        if int(eeg_norm) >= 1:                       # MATLAB: pwr in dB
            feats['frequency']['pwr'] = 10.0 * np.log10(np.maximum(PWR, 1e-20))
        feats['frequency'].update({k: np.round(v, 3)
                                   for k, v in band_pwr.items()})
        feats['qeeg'] = qeeg

        # ---- IAF (restingIAF) ------------------------------------------
        print('Attempting to find the individual alpha frequency (IAF) for '
              'each EEG channel...')
        try:
            pSum, pChans = resting_iaf(signals, [1, 30], [7, 14], 11, 5,
                                       fs=fs, cmin=1)
            feats['frequency']['IAF_mean'] = round(float(pSum['cog']), 3)
            feats['frequency']['IAF'] = np.round(pChans['cogs'], 3)
            if np.isfinite(pSum['cog']):
                print(f'Mean IAF across all channels: {pSum["cog"]:g}')
        except Exception as exc:
            print(f'IAF estimation failed: {exc}')
            feats['frequency']['IAF'] = np.full(n_chan, np.nan)

        # ---- Alpha asymmetry -------------------------------------------
        if n_chan > 1:
            alpha_pwr = PWR[:, (f >= bands[2, 0]) & (f <= bands[2, 1])].mean(
                axis=1, where=np.isfinite(PWR[:, (f >= bands[2, 0]) &
                                              (f <= bands[2, 1])]))
            alpha_pwr = np.nanmean(
                PWR[:, (f >= bands[2, 0]) & (f <= bands[2, 1])], axis=1)
            tot_pwr = np.nanmean(PWR, axis=1)
            asy, pairLabels, pairNums = compute_asymmetry(
                alpha_pwr, asy_norm, params['chanlocs'], False,
                tot_pwr if asy_norm else None)
            feats['frequency']['asymmetry'] = np.round(asy, 3)
            feats['frequency']['asymmetry_pairs_labels'] = list(pairLabels)
            feats['frequency']['asymmetry_pairs_num'] = pairNums
        else:
            print('WARNING: only one EEG channel: alpha asymmetry skipped.')

    # ---- Nonlinear domain -------------------------------------------------
    if params.get('eeg_nonlinear'):
        print('Computing EEG nonlinear features...')
        if fs > 100:
            new_fs = 90
            print(f'Resampling EEG data to {new_fs} Hz...')
            signals_use, fs_use = _resample_90(signals, fs)
        else:
            signals_use, fs_use = signals, fs
            new_fs = fs_use
        FD = np.empty(signals_use.shape[0])
        FE = np.empty(signals_use.shape[0])
        for iC in range(signals_use.shape[0]):
            print(f'  channel {iC + 1}...')
            fd, _sd = fractal_volatility(signals_use[iC])
            fe, _p = compute_fe(signals_use[iC], 2, 0.15, 2, 1)
            FD[iC], FE[iC] = fd, fe
        feats['nonlinear'] = {'FD': FD, 'FE': FE}

    return feats, params