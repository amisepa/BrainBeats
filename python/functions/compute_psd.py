"""compute_psd - Welch power spectral density of each channel.

Python port of BrainBeats functions/compute_psd.m: scipy.signal.welch with
MATLAB-pwelch semantics (hamming/hann/blackman/rectwin taper, overlap in %,
nfft = next power of 2 of the window length, one-sided 'density' PSD;
MATLAB's overlap in samples = detrended windows scaled by scipy's 'nperseg
/ noverlap' when overlap < 50... both handled below). Returns
(pwr, pwr_db, f) with pwr (channels x freqs) inside fRange.

Copyright (C) Cedric Cannard, 2021 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import get_window, welch

__all__ = ['compute_psd', 'nextpow2']


def nextpow2(n):
    """MATLAB nextpow2: smallest p with 2**p >= |n|."""
    n = abs(float(n))
    p = 0
    while 2 ** p < n:
        p += 1
    return p


def compute_psd(eegData, winSize=None, taperM='hamming', overlap=50,
                nfft=None, Fs=None, fRange=None, type='psd', useGPU=False):
    """eegData: (channels x samples). winSize: window length in SAMPLES
    (MATLAB convention: the callers pass EEG.srate*winlen_s). taperM:
    'hamming'|'hann'|'blackman'|'rectwin'. overlap: percent (default 50).
    nfft: default next power of 2 of winSize. fRange: [low high] Hz, default
    the full Nyquist band scaled like MATLAB ([1/nyquist, nyquist]).
    type: 'psd' (uV^2/Hz) or 'power' (PSD x window ENB: uV^2).
    """
    if Fs is None:
        raise ValueError('compute_psd: you need to provide the sampling rate Fs')
    data = np.asarray(eegData, float)
    one_d = data.ndim == 1
    if one_d:
        data = data[None, :]
    if not winSize:
        winSize = int(Fs * 2)
    if taperM is None:
        taperM = 'hamming'
    if not overlap:
        overlap = 50
    ov_samples = int(winSize / (100.0 / overlap))
    if fRange is None:
        nyq = Fs / 2.0
        fRange = [1.0 / nyq, nyq]
    if not nfft:
        nfft = 2 ** nextpow2(winSize)
    taper_names = {'hamming': 'hamming', 'hann': 'hann', 'blackman': 'blackman',
                   'rectwin': ('rect',)}
    win_name = taper_names[str(taperM).lower()]
    scale_map = {'psd': 'density', 'power': 'spectrum'}
    pwr_list, f = None, None
    for iChan in range(data.shape[0]):
        f, pxx = welch(data[iChan], fs=Fs, window=win_name,
                       nperseg=winSize, noverlap=ov_samples, nfft=nfft,
                       detrend='constant', scaling=scale_map[str(type).lower()],
                       return_onesided=True)
        if pwr_list is None:
            pwr_list = np.empty((data.shape[0], pxx.size))
        pwr_list[iChan] = pxx
    if str(type).lower() == 'power':
        # MATLAB pwelch 'power' option returns the PSD scaled by the ENB of
        # the window; scipy 'spectrum' scaling divides by the ENB instead:
        # correct it back so band integrals match the density integral.
        from scipy.signal.windows import get_window
        w = get_window(win_name, winSize)
        enb = winSize * (w ** 2).sum() / (w.sum() ** 2)
        pwr_list *= enb
    pwr = pwr_list
    keep = (f >= fRange[0]) & (f <= fRange[1])
    f = f[keep]
    pwr = pwr[:, keep]
    pwr_db = 10.0 * np.log10(np.maximum(pwr, np.finfo(float).tiny))
    if one_d:
        pwr, pwr_db = pwr[0], pwr_db[0]
    return pwr, pwr_db, f