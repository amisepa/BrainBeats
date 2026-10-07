"""run_features - HRV + EEG features mode driver (BrainBeats mode 2/3).

Mirrors the brainbeats_process.m features branch: given raw signals (ECG
for HRV, EEG for the EEG features), clean the NN series and compute the
HRV/EEG features, then optionally gather them into tables. Also exposes
the pop_ function expected by the eegprep-style plugin shell.

Copyright (C) Cedric Cannard, 2023 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import numpy as np

from .get_eeg_features import get_eeg_features
from .get_hrv_features import get_hrv_features

__all__ = ['run_features', 'pop_run_features']


def _clean_rr_series(rr, art_thresh=1.25, ectopic_thresh=0.85):
    """Drop RR intervals outside [1-art-thresh, 1+art-thresh] x median.
    (Same filter as brainbeats_process.m before the HRV features.)"""
    rr = np.asarray(rr, float).ravel()
    valid = np.isfinite(rr) & (rr > 0.3) & (rr < 3.0)
    rr = rr[valid]
    med = float(np.median(rr))
    ok = (rr > med * (1 - art_thresh)) & (rr < med * (1 + art_thresh))
    return rr[ok]


def run_features(ECG=None, EEG=None, NN=None, NN_times=None, params=None):
    """Compute HRV and/or EEG features.

    ECG: (1 x n) ECG signal uV; NN/NN_times override ECG (pre-detected RR).
    EEG: (nEEG x n) data uV; params: dict with keys:
      fs (REQUIRED), chanlocs (for EEG asymmetry), hrv_time, hrv_frequency,
      hrv_nonlinear, hrv_spec ('LombScargle_norm'/'LombScargle'/'welch'/
      'fft'), hrv_norm, hrv_overlap, eeg_time, eeg_frequency, eeg_nonlinear,
      eeg_frange, eeg_wintype, eeg_winlen, eeg_winoverlap, eeg_freqbounds
      ('conventional'/'individualized'), eeg_norm, asy_norm.
    Returns (Features, params): Features = {'HRV':..., 'EEG':...} with only
    the requested domains present.
    """
    params = dict(params) if params else {}
    if ECG is None and NN is None and EEG is None:
        raise ValueError("run_features: provide ECG/NN and/or EEG signals.")

    Features = {}
    if (params.get('hrv_time') or params.get('hrv_frequency')
            or params.get('hrv_nonlinear')):
        fs = params.get('fs') or params.get('ecg_fs')
        if NN is None:
            if ECG is None:
                raise ValueError('run_features: HRV needs ECG or NN input.')
            from ._sources import load_module            # validated assets
            get_rr = load_module('get_rr_v25.py').get_rr
            if not fs:
                fs = 250
            r_peaks, info = get_rr(np.asarray(ECG, float).ravel(), fs,
                                   drop_first=False)
            NN = np.diff(r_peaks) / fs
            NN_times = r_peaks[1:] / fs
        NN = np.asarray(NN, float).ravel()
        NN_times = np.asarray(NN_times, float).ravel()
        NN, NN_times = _clean_rr_series_pair(NN, NN_times)
        HRV, params = get_hrv_features(NN, NN_times, params)
        Features['HRV'] = HRV

    if params.get('eeg_time') or params.get('eeg_frequency') \
            or params.get('eeg_nonlinear'):
        if EEG is None:
            raise ValueError('run_features: EEG features need EEG input.')
        if not params.get('fs'):
            raise ValueError("run_features: params['fs'] is required.")
        EEGf, params = get_eeg_features(np.asarray(EEG, float), params)
        Features['EEG'] = EEGf

    return Features, params


def _clean_rr_series_pair(NN, NN_times):
    """Keep NN in [0.3, 3] s and within 25% of the median (clean_rr parity)."""
    ok = np.isfinite(NN) & (NN > 0.3) & (NN < 3.0)
    NN, NN_times = NN[ok], NN_times[ok]
    med = float(np.median(NN)) if NN.size else np.nan
    ok = np.isfinite(NN) & (NN > med * 0.75) & (NN < med * 1.25)
    return NN[ok], NN_times[ok]


def pop_run_features(EEG, *args, **kwargs):
    """eegprep-style entry: compute features on EEG.data / ECG channel and
    store them in EEG.brainbeats.features."""
    opts = kwargs.get('options') or kwargs
    data = EEG['data']
    fs = float(EEG['srate'])
    params = {**opts, 'fs': fs} if isinstance(opts, dict) else {'fs': fs}
    ECG = None
    chanlocs = EEG.get('chanlocs')
    params.setdefault('chanlocs', chanlocs)
    # ECG channel of choice, if any
    labels = [str(c.get('labels', '')) if isinstance(c, dict) else ''
              for c in (chanlocs or [])]
    for want in ('ECG', 'EKG'):
        if want in [lab.upper() for lab in labels]:
            ECG = data[labels.index(want if want in labels
                                    else [l for l in labels
                                          if l.upper() == want][0])]
            break
    params.setdefault('hrv_time', True)
    params.setdefault('hrv_frequency', True)
    params.setdefault('hrv_nonlinear', True)
    params.setdefault('eeg_time', True)
    params.setdefault('eeg_frequency', True)
    params.setdefault('eeg_nonlinear', True)
    params['eeg'] = params.get('eeg', True) and True
    Features, params = run_features(ECG=ECG, EEG=data, params=params)
    EEG.setdefault('brainbeats', {})['features'] = Features
    print('run_features OK: %d domain(s) computed.' % len(Features))
    return EEG, Features