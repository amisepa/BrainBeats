"""Compare our PPG beat train with the reference NN series by TIME (not a
greedy interval walk): ref beat times = start + cumsum(NN); ours = peaks/fs.
Find the offset that maximises matched beats (within 2 samples) and report
how many of the 301 reference beats our train reproduces.
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

HEP = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
ref_NN = np.asarray(HEP['brainbeats']['preprocessings']['NN'], float)
n_ref = ref_NN.size                 # 301

print('ref NN[:5] (s):', ref_NN[:5], 'median', np.median(ref_NN))
our_train = np.load(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/tests/our_train_ppg.npy')
print('our PPG kept beats:', our_train.size, 'NN[:5] (s):', np.diff(our_train/250.0)[:5])
our_t = (our_train - 1) / 250.0                     # 0-based s
ref_t = np.cumsum(ref_NN)                           # intervals -> beat times (s)

# offset scan: ref_t + off == our_t
best = (0, -1)
for off in np.arange(-1.0, 1.0, 0.001):
    diffs = np.abs(our_t[None, :] - (ref_t[:, None] + off))
    match = (diffs.min(axis=1) < 0.0081).sum()
    if match > best[1]:
        best = (off, match)
off, match = best
print(f'offset {off:.3f} s matches {match}/{n_ref} reference beats (<8.1 ms)')
# ECG train for comparison
fs = 250.0
ECG = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset.set')
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(ECG['chanlocs']).ravel()]
ecg = np.asarray(ECG['data'])[labels.index('ECG')].astype(float)
peaks_full, info = get_rr(ecg, fs, drop_first=False)
ecg_t = (peaks_full[1:] - 1) / 250.0
best2 = (0, -1)
for off2 in np.arange(-1.0, 1.0, 0.001):
    diffs = np.abs(ecg_t[None, :] - (ref_t[:, None] + off2))
    match2 = (diffs.min(axis=1) < 0.0081).sum()
    if match2 > best2[1]:
        best2 = (off2, match2)
print('ECG train matches ref PPG-derived NN beats:', best2[1], 'of', n_ref,
      f'at offset {best2[0]:.3f} s (pulse transit ~{best2[0]*1000:.0f} ms)')
# also: ECG vs our PPG train offset (expected PAT ~ tens of ms)
best3 = (0, -1)
for off3 in np.arange(-0.5, 0.5, 0.001):
    diffs = np.abs(ecg_t[None, :] - (our_t[:, None] + off3))
    m3 = (diffs.min(axis=1) < 0.0081).sum()
    if m3 > best3[1]:
        best3 = (off3, m3)
print('ECG vs our-PPG train:', best3[1], 'matches at offset', best3[0])