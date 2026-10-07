"""TF-core sanity: run compute_hep_tf and compute the HEP by hand on the
MATLAB-saved EPOCHS of dataset_HEP.set. Zero alignment uncertainty: the beats
are the epoch triggers, the data is the MATLAB-cleaned signal.

Also cross-checks the grid: reference hrsp.times vs our tIdx grid, and the
reference hep_times (224 pts) vs our hepIdx (225 pts at these params).
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
import eegprep
from functions.compute_hep_tf import compute_hep_tf

HEP = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
hrsp_ref = HEP['brainbeats']['hrsp']
data = np.asarray(HEP['data'], float)                # (63, 224, 267)
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(HEP['chanlocs']).ravel()]
fs = HEP['srate']
roi = ['F1','F2','F3','F4','Fz','FC1','FC2','FC3','FC4','FCz',
       'C1','C2','C3','C4','Cz']
ridx = [i for i, l in enumerate(labels) if l in roi]
print('n roi chans:', len(ridx))

# 1) ROI HEP by hand: mean over beats then mean over ROI channels
hep_hand = data[ridx].mean(axis=0).mean(axis=1)      # (224,)
hep_times_ref = np.asarray(hrsp_ref['hep_times'], float)
print('ref hep_times grid:', hep_times_ref[0], '->', hep_times_ref[-1],
      'n', hep_times_ref.size, 'step', np.unique(np.diff(hep_times_ref)))
print('epo times:', np.asarray(HEP['times']).ravel()[0],
      np.asarray(HEP['times']).ravel()[-1])

# 2) our TF core on the EPOCH-MEAN signal directly (a single channel spanning
#    224 samples, fs 250): beats all at sample 0 of this concatenated view
sig_roi = hep_hand
# build one long continuous signal from the 267 cleaned epochs to avoid edge
# effects in the wavelet convolution: concat with no overlap
longsig = data[ridx].mean(axis=0).T.ravel()          # (224*267,)
beats = (np.arange(data.shape[2]) * 224 + 1).astype(np.int64)  # 1-based
tf, _ = compute_hep_tf(longsig[None, :], fs, beats,
                       (hep_times_ref[0], hep_times_ref[-1] + 1),
                       hep_times=hep_times_ref, n_surr=0)
hep_py = tf['hep'][0]
r = np.corrcoef(hep_py, hep_hand)[0, 1]
d = np.abs(hep_py - hep_hand)
print(f'TF-core HEP vs MATLAB epoch HEP: corr {r:.6f}  max|diff| {d.max():.3e} uV')

# 3) HRSP core on the LONG concatenated signal, all channels, our grids
# (this bypasses cleaning: the stored epochs ARE the cleaned data)
long_all = data.reshape(data.shape[0], -1)
beats_all = (np.arange(data.shape[2]) * 224 + 1).astype(np.int64)
win = (-300.0, 600.0)
times_ref = np.asarray(hrsp_ref['times'], float)
freqs_ref = np.asarray(hrsp_ref['freqs'], float)
print('ref freqs:', freqs_ref.size, freqs_ref[:6], '...', freqs_ref[-3:])
tf_all, _ = compute_hep_tf(long_all, fs, beats_all, win, hep_times=hep_times_ref,
                           freqs=freqs_ref, n_surr=0)
hrsp_py = tf_all['hrsp']
hrsp_ref_arr = np.asarray(hrsp_ref['hrsp'], float)
print('shapes:', hrsp_py.shape, hrsp_ref_arr.shape)
r_b = np.corrcoef(hrsp_py.ravel(), hrsp_ref_arr.ravel())[0, 1]
d_b = np.abs(hrsp_py - hrsp_ref_arr)
print(f'TF-core HRSP vs reference: corr {r_b:.6f} max|dB| {d_b.max():.3f} '
      f'median {np.median(d_b):.4f}')
hepc_py = tf_all['hrpc']
hepc_ref = np.asarray(hrsp_ref['hepc'], float)
r_c = np.corrcoef(hepc_py.ravel(), hepc_ref.ravel())[0, 1]
print(f'TF-core HRPC vs reference: corr {r_c:.6f} max {np.abs(hepc_py - hepc_ref).max():.4f}')