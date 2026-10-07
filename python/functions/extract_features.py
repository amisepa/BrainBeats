"""extract_features - Gather HRV and EEG features into one-row tables.

Python port of BrainBeats functions/extract_features.m: takes the
brainbeats_process features dict (keys 'HRV' and/or 'EEG') and returns
(hrv_df, eeg_df) pandas DataFrames -- one row, one column per feature
(channel features flattened to one column per channel in channel order).
Empty (but non-None) DataFrames when a domain is absent.

Copyright (C) Cedric Cannard, 2023 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ['extract_features']


def _flat(x):
    x = np.asarray(x, float).ravel()
    return x


def _add_var(cols, x, name):
    """Append x as one column per element (MATLAB add_var subfunction)."""
    vals = _flat(x)
    if vals.size == 1:
        cols[f'{name}'] = [vals[0]]
    else:
        for i, v in enumerate(vals):
            cols[f'{name}_{i + 1}'] = [v]


def extract_features(Features):
    hrv_cols, eeg_cols = {}, {}

    HRV = Features.get('HRV') if Features else None
    if HRV:
        # Time
        tm = HRV.get('time')
        if tm:
            for k, v in tm.items():
                hrv_cols[f'{k}-HRV'] = [float(v)]
        # Frequency (mean across windows)
        freq = HRV.get('frequency')
        if freq:
            for band, name in (('ulf', 'ULF-HRV'), ('vlf', 'VLF-HRV'),
                               ('lf', 'LF-HRV'), ('hf', 'HF-HRV'),
                               ('lfhf', 'LF/HF')):
                if band in freq:
                    hrv_cols[name] = [float(np.mean(freq[band]))]
        # Nonlinear
        nl = HRV.get('nonlinear')
        if nl:
            poc = nl.get('Poincare')
            if poc:
                hrv_cols['Poincaré: SD1'] = [float(poc['SD1'])]
                hrv_cols['Poincaré: SD2'] = [float(poc['SD2'])]
                hrv_cols['SD1/SD2'] = [float(poc['SD1SD2'])]
            if 'PRSA_AC' in nl and 'PRSA_DC' in nl:
                hrv_cols['PRSA_AC'] = [float(nl['PRSA_AC'])]
                hrv_cols['PRSA_DC'] = [float(nl['PRSA_DC'])]
            if 'FE' in nl:
                hrv_cols['HRV-FE'] = [float(nl['FE'])]
            if 'FD' in nl:
                hrv_cols['HRV-FD'] = [float(nl['FD'])]
            if 'MFE' in nl:                       # multiscale FE (if present)
                mfe = np.asarray(nl['MFE'], float)
                hrv_cols['HRV-MFE_peak'] = [int(np.argmax(mfe)) + 1]
                hrv_cols['HRV-MFE_auc'] = [float(np.trapezoid(mfe))]

    EEG = Features.get('EEG') if Features else None
    if EEG:
        tm = EEG.get('time')
        if tm:
            for k, v in tm.items():
                _add_var(eeg_cols, v, f'EEG-{k}')
        freq = EEG.get('frequency')
        if freq:
            for band in ('delta', 'theta', 'alpha', 'beta', 'gamma'):
                if band in freq:
                    _add_var(eeg_cols, freq[band], f'EEG-{band}')
            if 'IAF' in freq:
                _add_var(eeg_cols, freq['IAF'], 'EEG-IAF')
            if 'IAF_mean' in freq:
                _add_var(eeg_cols, freq['IAF_mean'], 'EEG-IAF_mean')
            asy = freq.get('asymmetry')
            pairs = freq.get('asymmetry_pairs_labels')
            if asy is not None and pairs is not None:
                for i, (a, p) in enumerate(zip(np.asarray(asy).ravel(), pairs)):
                    eeg_cols[f'Asy ({p})'] = [float(a)]
        qe = EEG.get('qeeg')
        if qe:
            for k, v in qe.items():
                _add_var(eeg_cols, v, f'qEEG-{k}')
        nl = EEG.get('nonlinear')
        if nl:
            for k, v in nl.items():
                _add_var(eeg_cols, v, f'EEG-{k}')

    return pd.DataFrame(hrv_cols), pd.DataFrame(eeg_cols)