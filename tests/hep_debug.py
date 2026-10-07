"""Isolate the HEP mismatch: run our TF core's HEP path AND a hand-rolled
epoch mean on the MATLAB-cleaned stored epochs of dataset_HEP.set. If the
hand mean matches hep_hand (corr 1.0) but the TF core doesn't, the bug is in
_mean_epochs/valid_beats; if both mismatch, the mapping (concat trick /
win/ hep_times) is wrong.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
import eegprep
import functions.compute_hep_tf as C

HEP = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
hrsp_ref = HEP['brainbeats']['hrsp']
data = np.asarray(HEP['data'], float)                 # (63, 224, 267)
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(HEP['chanlocs']).ravel()]
roi = ['F1','F2','F3','F4','Fz','FC1','FC2','FC3','FC4','FCz',
       'C1','C2','C3','C4','Cz']
ridx = [i for i, l in enumerate(labels) if l in roi]
hep_hand = data[ridx].mean(axis=0).mean(axis=1)       # (224,)
hep_times_ref = np.asarray(hrsp_ref['hep_times'], float)
longsig = data[ridx].mean(axis=0).T.ravel()           # (224*267,) = epochs concatenated
beats = (np.arange(data.shape[2]) * 224 + 1).astype(np.int64)  # 1-based

tf, _ = C.compute_hep_tf(longsig[None, :], 250.0, beats,
                         (hep_times_ref[0], hep_times_ref[-1] + 1),
                         hep_times=hep_times_ref, n_surr=0, tf=False)
hep_py = tf['hep'][0]
print('TF core kept', tf['nBeats'], 'beats (of', beats.size, ')')
r_py = np.corrcoef(hep_py, hep_hand)[0, 1]
print(f'TF-core HEP vs stored-truth:  corr {r_py:.6f}')

# hand-rolled mean with explicit 1-based -> 0-based mapping: sig[(b-1) + idx]
idx0 = np.rint(hep_times_ref / 1000.0 * 250.0).astype(int)   # -75..148
man = np.stack([longsig[(b - 1) + idx0] for b in beats]).mean(axis=0)
r_man = np.corrcoef(man, hep_hand)[0, 1]
d_man = np.abs(man - hep_hand)
print(f'hand mean vs stored-truth:    corr {r_man:.8f}  max|diff| {d_man.max():.3e}')
# where do they differ (if at all)?
print('hep_hand[:4]', hep_hand[:4])
print('man      [:4]', man[:4])
print('hep_py   [:4]', hep_py[:4])
print('hep_py[-4:]', hep_py[-4:])
print('std hep_py vs man:', hep_py.std(), man.std())