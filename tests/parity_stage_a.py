"""Parity driver: reproduce the dataset_HEP.set analysis in Python.

Chain (matching the MATLAB reference run):
  pop_loadset -> pop_select(nochannel, PPG ECG) -> pop_eegfiltnew(0.5, 30)
  -> apply_car (full-rank CAR) -> clean_channels(corr .65, line 15 sd, win 5,
  maxBad .33) -> interpolate TP9 back (spherical) -> get the R-peak train via
  Pan-Tompkins (get_rr) + clean_rr -> beat rejections (window + IBI Grubbs)
  -> epochs at -950..+1250 ms -> find_badTrials -> Picard ICA at rank
  -> ICLabel -> icflag -> subcomp (heart .75) -> crop -300..600 ms
  -> T estimate -> clean continuous -> compute_hep_tf on the ROI mean and all
  channels (4:30 Hz, cycles 5, tstep 10 ms -> 12 ms effective at 250 Hz).
Compare: HRSP/HRPC/power vs reference_hep_tf.npz.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
import eegprep

BB = r'C:/Users/ccann/Documents/MATLAB/BrainBeats'
SAMPLE = os.path.join(BB, 'sample_data', 'dataset.set')

EEG = eegprep.pop_loadset(SAMPLE)
EEG = eegprep.pop_select(EEG, 'nochannel', ['PPG', 'ECG'])
EEG = eegprep.pop_eegfiltnew(EEG, 0.5, 30)

# full-rank CAR (apply_car.m, sans the reference recovery)
X = np.asarray(EEG['data'], float)
rank_b = int((np.linalg.eigvalsh(np.cov(X)) > 1e-7).sum())
X0 = np.vstack([X, np.zeros((1, X.shape[1]))])
X0 = X0 - X0.mean(axis=0, keepdims=True)
EEG['data'] = X0[:-1]
rank_a = int((np.linalg.eigvalsh(np.cov(EEG['data'])) > 1e-7).sum())
print('apply_car: rank before/after:', rank_b, rank_a)

# bad channels (clean_rawdata port in eegprep)
out = eegprep.clean_channels(dict(EEG), corr_threshold=0.65, noise_threshold=15,
                             window_len=5, max_broken_time=0.33)
badEEG = out[0] if isinstance(out, tuple) else out
mask = np.asarray(badEEG['etc']['clean_channel_mask'], bool)
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(EEG['chanlocs']).ravel()]
bad_ch = [labels[i] for i in np.flatnonzero(~mask)]
print('bad channels:', bad_ch)
EEG = badEEG
# interpolate back (clean_eeg: removed channels re-interpolated to orichanlocs)
bad_idx = np.flatnonzero(~mask) + 1          # 1-based for pop_interp
out = eegprep.pop_interp(EEG, bad_idx.tolist(), 'spherical')
EEG = out[0] if isinstance(out, tuple) else out
print('interpolated; nbchan', EEG['nbchan'])

np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'parity_data.npz'),
         data=np.asarray(EEG['data']),
         lat=[ev['latency'] for ev in eegprep.pop_loadset(
             os.path.join(BB, 'sample_data', 'dataset_HEP.set'))['event']])
print('stage-A done: cleaned continuous data + 267 reference R-peak latencies saved.')