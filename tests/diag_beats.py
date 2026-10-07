"""Diagnose the beat-count difference on the sample dataset:
MATLAB NN has 301 intervals; our get_rr + clean_rr gave 305 clean peaks.
Steps: count raw detections, then clean_rr with interpolate_missing=True (the
MATLAB brainbeats_process call), report removed/synthetic/kept.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/HEP_neurofeedback/python')
import eegprep

EEG = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset.set')
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(EEG['chanlocs']).ravel()]
ecg = np.asarray(EEG['data'])[labels.index('ECG')].astype(float)
fs = EEG['srate']

from get_rr_v25 import get_rr
from rr_cleaning_v25 import clean_rr, qrs_bandpass

peaks_full, info = get_rr(ecg, fs, drop_first=False)
print('raw detections:', peaks_full.size)
filt = qrs_bandpass(ecg, fs)
peaks_oneb = peaks_full[1:]
rr = np.diff(peaks_full) / fs
peak_amp = filt[peaks_oneb - 1]
nn, npk, idx_bad, info_rr = clean_rr(
    rr, peaks_oneb, fs, peak_amp=peak_amp, ecg_signal=ecg,
    sig_t=np.arange(ecg.size) / fs, interpolate_missing=True, verbose=True)
print('clean_rr: nn', nn.size, 'npk', npk.size,
      'synthetic', int(np.sum(~np.isfinite(npk))),
      'idx_bad', int(idx_bad.sum()))
kept_real = npk[np.isfinite(npk) & ~idx_bad]
print('kept real (excl synthetic):', int(kept_real.size))
print('MATLAB reference: NN=301')
# first/last beat times
print('first kept beat:', kept_real[0], 'last:', kept_real[-1])
print('raw first beat (1-based):', peaks_full[0], 'last:', peaks_full[-1])