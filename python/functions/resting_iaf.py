"""resting_iaf - Individual alpha frequency (IAF; PAF and CoG) of resting EEG.

Python port of the restingIAF package (Corcoran et al. 2018,
github.com/corcorana/restingIAF) as called by BrainBeats
get_eeg_features.m: restingIAF(signals, nChan, 1, [1 30], fs, [7 14], 11, 5)
-> per-channel alpha centre of gravity (CoG) and its cross-channel mean.

Also exposes get_freqBounds (BrainBeats functions/get_freqBounds.m, the
restingIAF method adapted to any band).

Copyright (C) Andrew W. Corcoran 2016-2018 (MATLAB);
BrainBeats wrapper Cedric Cannard 2023; Python port 2026-10.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import savgol_filter, welch
from scipy.signal.windows import hamming

from .compute_psd import nextpow2

__all__ = ['resting_iaf', 'get_freqBounds', 'findF1', 'findF2']


def _nearest(f, value):
    """dsearchn(f, value): index of the nearest frequency bin."""
    return int(np.argmin(np.abs(np.asarray(f) - value)))


def _pwelch_matlab(x, Fs, taper='hamming', tlen=None, tover=None, nfft=None):
    """MATLAB pwelch(x, window, noverlap, nfft, fs): one-sided density PSD.

    Defaults: noverlap = 50% of window; nfft = nextpow2 window; hamming.
    """
    x = np.asarray(x, float).ravel()
    if tlen is None:
        tlen = int(Fs * 4)
    tlen = int(min(tlen, x.size))
    if tover is None:
        tover = tlen // 2
    if not nfft:
        nfft = 2 ** nextpow2(tlen)
    win = hamming(tlen) if str(taper) == 'hamming' else None
    f, pxx = welch(x, fs=Fs, window=win, nperseg=tlen, noverlap=tover,
                   nfft=nfft, detrend='constant', scaling='density',
                   return_onesided=True)
    return f, pxx


def _sgf_diff(x, Fw, poly, Fs, tlen):
    """restingIAF sgfDiff: scipy savgol_filter deriv mirrors MATLAB
    sgolay + conv, for derivative orders 0, 1, 2."""
    dt = Fs / float(tlen)                    # bin width of the pwelch grid
    d = [savgol_filter(x, int(Fw), int(poly), deriv=p, delta=dt)
         for p in (0, 1, 2)]
    # MATLAB's fact(p)/(-dt)^p * g(:,p+1) + conv reduces to scipy's
    # deriv=p/delta=dt (the (-dt)^p sign compensates conv's kernel reversal),
    # i.e. d1 = +dPSD/dHz ascending into the peak.
    return d[0], d[1], d[2]


def _lower_upper_zero_crossings(d0, d1, f, lo, hi):
    """Downward crossings of d1 in [lo-1, hi+1] (restingIAF peakBounds)."""
    rows = []
    for k in range(max(lo - 1, 0), min(hi + 1, len(d1) - 1)):
        if np.sign(d1[k]) > np.sign(d1[k + 1]):       # downward crossing
            maxim = k if d0[k] >= d0[k + 1] else k + 1
            rows.append([len(rows) + 1, maxim, float(f[maxim]), d0[maxim]])
    return rows


def _peak_selection(rows, minPow, mdiff):
    """restingIAF: pick (peakBin, peakF, subBin) from candidate crossings."""
    if not rows:
        return None, None, None
    if len(rows) == 1:
        r = rows[0]
        if np.log10(r[3]) > minPow[r[1]]:
            return r[1], r[2], None
        return None, None, None
    rows = sorted(rows, key=lambda r: -r[3])
    r0, r1 = rows[0], rows[1]
    if np.log10(r0[3]) > minPow[r0[1]] and r0[3] * (1 - mdiff) > r1[3]:
        return r0[1], r0[2], None
    if np.log10(r0[3]) > minPow[r0[1]]:
        return None, None, r0[1]
    return None, None, None


def findF1(f, d0, d1, negZ, minPow, slen, bin_):
    """Lower alpha bound (restingIAF findF1 / BrainBeats findF1.m).""" ""
    negZ = [list(r) for r in negZ]
    if len(negZ) > 1:
        negZ = sorted(negZ, key=lambda r: r[2])
        leftPeak = bin_
        for z in range(len(negZ)):
            if (np.log10(negZ[z][3]) > minPow[negZ[0][1]]
                    or negZ[z][3] > 0.5 * d0[bin_]):
                leftPeak = negZ[z][1]
                break
    else:
        leftPeak = bin_
    posZ1 = []
    for k in range(1, leftPeak):            # MATLAB 2:leftPeak-1 -> py 1..
        if k + 1 >= len(d1):
            break
        if np.sign(d1[k]) < np.sign(d1[k + 1]):
            trio = np.abs([d0[k - 1], d0[k], d0[k + 1]])
            minim = k - 1 + int(np.argmin(trio))
            posZ1.append([len(posZ1) + 1, minim, f[minim]])
        elif abs(d1[k]) < 1 and np.all(np.abs(d1[k + 1:min(k + slen, len(d1))]) < 1):
            posZ1.append([len(posZ1) + 1, k, f[k]])
    if not posZ1:
        return np.nan, np.nan
    if len(posZ1) == 1:
        return posZ1[0][1], posZ1[0][2]
    posZ1 = sorted(posZ1, key=lambda r: -r[2])
    return posZ1[0][1], posZ1[0][2]


def findF2(f, d0, d1, negZ, minPow, slen, bin_):
    """Upper alpha bound (restingIAF findF2)."""
    negZ = [list(r) for r in negZ]
    if len(negZ) > 1:
        negZ = sorted(negZ, key=lambda r: -r[2])
        rightPeak = bin_
        for z in range(len(negZ)):
            if (np.log10(negZ[z][3]) > minPow[negZ[0][1]]
                    or negZ[z][3] > 0.5 * d0[bin_]):
                rightPeak = negZ[z][1]
                break
    else:
        rightPeak = bin_
    posZ2 = []
    for k in range(rightPeak + 1, len(d1) - slen):
        if np.sign(d1[k]) < np.sign(d1[k + 1]):
            trio = np.abs([d0[k - 1], d0[k], d0[k + 1]])
            minim = k - 1 + int(np.argmin(trio))
            posZ2.append([len(posZ2) + 1, minim, f[minim]])
        elif abs(d1[k]) < 1 and np.all(np.abs(d1[k + 1:k + slen]) < 1):
            posZ2.append([len(posZ2) + 1, k, f[k]])
    if not posZ2:
        return np.nan, np.nan
    return posZ2[0][1], posZ2[0][2]         # first estimate only


def _peak_bounds(d0, d1, d2, f, w, minPow, mdiff, fres):
    """restingIAF peakBounds: (peakF, pos1, pos2, f1, f2, inf1, inf2, Q, Qf)"""
    lo = _nearest(f, w[0])
    hi = _nearest(f, w[1])
    rows = _lower_upper_zero_crossings(d0, d1, f, lo, hi)
    peakBin, peakF, subBin = _peak_selection(rows, minPow, mdiff)
    slen = int(round(1.0 / fres))
    if peakBin is None and subBin is None:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    start = peakBin if peakBin is not None else subBin
    f1, posZ1 = findF1(f, d0, d1, rows, minPow, slen, start)
    f2, posZ2 = findF2(f, d0, d1, rows, minPow, slen, start)
    if peakBin is None:
        return np.nan, posZ1, posZ2, f1, f2, np.nan, np.nan, np.nan, np.nan
    # inflection points around the peak (ascending edge = last downward
    # d2 crossing before the peak; descending = first upward crossing after)
    inf1, inf2 = np.nan, np.nan
    min1 = min2 = np.nan
    for k in range(0, peakBin - 1):
        if np.sign(d2[k]) > np.sign(d2[k + 1]):
            min1 = k if abs(d2[k]) <= abs(d2[k + 1]) else k + 1
    if np.isfinite(min1):
        inf1 = f[min1]
    for k in range(peakBin + 1, len(d2) - 1):
        if np.sign(d2[k]) < np.sign(d2[k + 1]):
            min2 = k if abs(d2[k]) <= abs(d2[k + 1]) else k + 1
            inf2 = f[min2]
            break
    if np.isfinite(min1) and np.isfinite(min2) and min2 > min1:
        Q = float(np.trapezoid(d0[min1:min2 + 1], f[min1:min2 + 1]))
        Qf = Q / (min2 - min1)
    else:
        Q = Qf = np.nan
    return peakF, posZ1, posZ2, f1, f2, inf1, inf2, Q, Qf


def resting_iaf(data, fRange=(1, 30), w=(7, 14), Fw=11, k=5, *, fs=None,
                mpow=1.0, mdiff=0.20, taper='hamming', tlen=None, tover=None,
                nfft=None, norm=True, cmin=1, window_size=None):
    """IAF of each channel (restingIAF semantics, BrainBeats call defaults).

    data: (nChan x samples). fs: sample rate (required). Returns
    (pSum, pChans) with pSum.cog / pSum.paf the cross-channel means (NaN
    below cmin) and pChans['gravs'] / ['peaks'] the per-channel CoG / PAF.
    """
    X = np.asarray(data, float)
    n_chan = X.shape[0]
    if fs is None:
        raise ValueError('resting_iaf: fs is required')
    if tlen is None:
        tlen = int(fs * 4)                       # MATLAB default tlen = Fs*4
    f_full, p_all = None, None
    per = []
    for iC in range(n_chan):
        row = X[iC]
        if np.isnan(row).any():
            per.append(None)                     # MATLAB: skip + trim
            continue
        f, pxx = _pwelch_matlab(row, fs, taper, tlen, tover, nfft)
        lo = _nearest(f, fRange[0])
        hi = _nearest(f, fRange[1])
        frex = np.arange(lo, hi + 1)             # dsearchn range, end-inclusive
        f = f[frex]
        pxx = pxx[frex]
        if norm:
            pxx = pxx / pxx.mean()
        # minPow: 1st-order log10 fit + mpow SD
        pfit = np.polyfit(f, np.log10(pxx), 1)
        yval = np.polyval(pfit, f)
        resid = np.log10(pxx) - yval
        sig = np.sqrt(np.sum(resid ** 2) / (len(f) - 2))    # polyfit sig (approx)
        del_ = sig * np.sqrt(1.0 / len(f) +
                             (f - f.mean()) ** 2 / np.sum((f - f.mean()) ** 2))
        minPow = yval + mpow * del_
        d0, d1, d2 = _sgf_diff(pxx, Fw, k, fs, tlen)
        exp_tlen = nextpow2(tlen)
        fres = fs / 2.0 ** exp_tlen
        out = _peak_bounds(d0, d1, d2, f, w, minPow, mdiff, fres)
        per.append({'f': f, 'pxx': pxx, 'minPow': minPow, 'd0': d0, 'd1': d1,
                    'd2': d2, 'iaw': None, 'gravs': None,
                    'peaks': out[0], 'f1': out[3], 'f2': out[4],
                    'Q': out[7], 'Qf': out[8]})
    f_full = next(p['f'] for p in per if p is not None)
    f1_f = np.array([p['f1'] if p else np.nan for p in per], float)
    f2_f = np.array([p['f2'] if p else np.nan for p in per], float)
    # f1/f2 are BIN indices: keep an integer copy for indexing (NaN where absent)
    f1 = np.where(np.isfinite(f1_f), np.round(f1_f).astype(int), -1)
    f2 = np.where(np.isfinite(f2_f), np.round(f2_f).astype(int), -1)
    d0s = np.column_stack([p['d0'] if p else np.full_like(f_full, np.nan)
                           for p in per])
    # chanGravs (trim off any NaNs = -1 sentinels); keep BIN indices ints
    trim_f1 = f1[f1 >= 0]
    trim_f2 = f2[f2 >= 0]
    iaw = (_nearest(f_full, float(f_full[trim_f1].mean())),
           _nearest(f_full, float(f_full[trim_f2].mean())))
    mf1, mf2 = iaw
    cogs = np.full(n_chan, np.nan)
    selG = np.isfinite(f1_f)                    # MATLAB sel = ~isnan(f1)
    for d in range(n_chan):
        if trim_f1.size == 0 or trim_f2.size == 0:
            break                                # all cogs NaN (MATLAB if/else)
        seg, fseg = d0s[mf1:mf2 + 1, d], f_full[mf1:mf2 + 1]
        cogs[d] = np.nansum(seg * fseg) / np.sum(seg)
    # chanMeans
    peaks = np.array([p['peaks'] if p else np.nan for p in per], float)
    qf = np.array([p['Qf'] if p else np.nan for p in per], float)
    selP = np.isfinite(peaks)
    with np.errstate(invalid='ignore'):
        wts = qf / np.nanmax(qf) if np.isfinite(qf).any() else qf
    if selP.sum() < cmin:
        paf = paf_std = mu_spec = np.nan
    else:
        wts = np.nan_to_num(wts)
        paf = (np.nansum(np.nan_to_num(peaks) * wts) / np.nansum(wts))
        paf_std = float(np.nanstd(peaks))
        mu_spec = np.nansum(np.nan_to_num(d0s) * wts[None, :], axis=1) / np.nansum(wts)
    pSum = {'paf': paf, 'pafStd': paf_std, 'muSpec': mu_spec,
            'cog': (float(np.nanmean(cogs)) if selG.sum() >= cmin else np.nan),
            'cogStd': (float(np.nanstd(cogs)) if selG.sum() >= cmin else np.nan),
            'pSel': int(selP.sum()), 'gSel': int(selG.sum()), 'iaw': iaw}
    pChans = {'f': f_full, 'cogs': cogs, 'peaks': peaks, 'qf': qf,
              'selG': selG, 'selP': selP, 'iaw': iaw,
              'd0': d0s}
    return pSum, pChans


# ---------------------------------------------------------------------------
# BrainBeats get_freqBounds.m (the method applied to any band)
# ---------------------------------------------------------------------------
def get_freqBounds(pwr, f, fs, w, winSize, mpow):
    """(bounds, peak) of a band from the spectrum shape: normalize, SGF
    differentiation, primary peak above the log10 1/f fit + mpow SD within
    w, bounds = minima either side. Returns ([lo hi], peak) in Hz (NaN where
    undetermined)."""
    f = np.asarray(f, float)
    p = np.asarray(pwr, float) / np.asarray(pwr, float).mean()
    fres = f[1] - f[0]
    Fw = max(int(2 * np.floor(2.69 / fres / 2) + 1), 5)
    k = min(5, Fw - 2)
    mdiff = 0.2
    pfit = np.polyfit(f, np.log10(p), 1)
    yval = np.polyval(pfit, f)
    resid = np.log10(p) - yval
    sig = np.sqrt(np.sum(resid ** 2) / (len(f) - 2))
    del_ = sig * np.sqrt(1.0 / len(f) + (f - f.mean()) ** 2 /
                         np.sum((f - f.mean()) ** 2))
    minPow = yval + mpow * del_
    d0, d1, d2 = _sgf_diff(p, Fw, k, fs, winSize)
    lo = _nearest(f, w[0])
    hi = _nearest(f, w[1])
    rows = []
    cnt = 0
    for kk in range(max(lo - 1, 0), min(hi + 1, len(d1) - 1)):
        if np.sign(d1[kk]) > np.sign(d1[kk + 1]):
            maxim = kk if d0[kk] >= d0[kk + 1] else kk + 1
            cnt += 1
            rows.append([cnt, maxim, f[maxim], d0[maxim]])
    if not rows:
        peak = np.nan
        rows = []
    elif len(rows) == 1:
        r = rows[0]
        peak = r[2] if np.log10(r[3]) > minPow[r[1]] else np.nan
    else:
        rows = sorted(rows, key=lambda r: -r[3])
        peak = (rows[0][2] if np.log10(rows[0][3]) > minPow[rows[0][1]]
                and rows[0][3] * (1 - mdiff) > rows[1][3] else np.nan)
    slen = int(round(1.0 / fres))
    if not rows:
        return [np.nan, np.nan], np.nan
    start = rows[0][1] if np.isfinite(peak) else (rows[0][1] if len(rows) > 1 else rows[0][1])
    f1, pos1 = findF1(f, d0, d1, rows, minPow, slen, start)
    f2, pos2 = findF2(f, d0, d1, rows, minPow, slen, start)
    return [pos1, pos2], peak