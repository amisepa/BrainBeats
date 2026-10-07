"""PPG beat train: our get_rr + clean_rr on the sample PPG channel.
The reference run used PPG (heart_channel_used = 'PPG', NN n=301).
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
sys.path.insert(0, r'C:/Users/ccann/Documents/HEP_neurofeedback/python')
import eegprep
from get_rr_v25 import get_rr
from rr_cleaning_v25 import clean_rr, qrs_bandpass

EEG = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset.set')
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(EEG['chanlocs']).ravel()]
ppg = np.asarray(EEG['data'])[labels.index('PPG')].astype(float)
fs = EEG['srate']
peaks_full, info = get_rr(ppg, fs, drop_first=False)
print('PPG raw detections:', peaks_full.size)
filt = qrs_bandpass(ppg, fs)
peaks_oneb = peaks_full[1:]
rr = np.diff(peaks_full) / fs
nn, npk, idx_bad, info_rr = clean_rr(
    rr, peaks_oneb, fs, peak_amp=filt[peaks_oneb - 1], ecg_signal=ppg,
    sig_t=np.arange(ppg.size) / fs, interpolate_missing=True)
npk = np.asarray(npk, float)
real = npk[np.isfinite(npk) & ~np.asarray(idx_bad, bool)]
print('PPG kept:', real.size, '(MATLAB reference: 301)')
nn_out = np.diff(real) / fs
ref_NN = None
HEP = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
ref_NN = np.asarray(HEP['brainbeats']['preprocessings']['NN'], float)
print('ref NN n', ref_NN.size, 'our NN n', nn_out.size)
# how many intervals match within 2 samples?
m = 0
i = 0
deltas = []
for iv_ref in ref_NN:
    while i < nn_out.size and nn_out[i] < iv_ref - 5.1e-3:
        i += 1
    if i < nn_out.size and abs(nn_out[i] - iv_ref) <= 5.1e-3:
        m += 1
        deltas.append(abs(nn_out[i] - iv_ref))
        i += 1
print('matched intervals:', m, 'of', ref_NN.size,
      'median delta', np.median(deltas) if deltas else 'NA')
np.save(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/tests/our_train_ppg.npy',
        real.astype(np.int64))