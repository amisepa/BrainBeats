"""Dump dataset_HEP.set preprocessings + the exact cleaning metadata, and test
whether a Python baseline_regression of the stored-epoch mean reproduces the
stored HEP (i.e. what the epochs were corrected with, if anything).
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
import eegprep

HEP = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
br = HEP['brainbeats']
pre = br['preprocessings']
print('preprocessings keys:', sorted(pre.keys()))
for k in sorted(pre.keys()):
    v = pre[k]
    if isinstance(v, dict):
        print(' ', k, '->', {kk: (np.asarray(vv).shape if hasattr(vv, 'shape') or
                                  isinstance(vv, (int, float, str, list)) else type(vv))
                             for kk, vv in v.items()})
    else:
        arr = np.asarray(v)
        print(' ', k, '->', arr.shape if arr.size > 1 else arr)

data = np.asarray(HEP['data'], float)
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(HEP['chanlocs']).ravel()]
roi = ['F1','F2','F3','F4','Fz','FC1','FC2','FC3','FC4','FCz',
       'C1','C2','C3','C4','Cz']
ridx = [i for i, l in enumerate(labels) if l in roi]
hep_hand = data[ridx].mean(axis=0).mean(axis=1)

# The reference .hep_times is 224 pts (-300..592). Is there any hep-like field
# to compare with? No roi struct, no hrsp.hep. So compare hep_hand with a
# regression variant of itself:
from functions.baseline_regression import baseline_regression
times_ms = np.asarray(hrsp_ref['hep_times'], float) if False else np.asarray(
    HEP['times'], float).ravel()
print('epoch times:', times_ms[0], '..', times_ms[-1])
try:
    ep_roi = data[ridx]
    out, beta, bl = baseline_regression(ep_roi, times_ms, (-150.0, -50.0))
    hep_br = out.mean(axis=0).mean(axis=1)
    r = np.corrcoef(hep_hand, hep_br)[0, 1]
    d = np.abs(hep_hand - hep_br).max()
    print(f'corr(stored, our baseline_regression of stored) {r:.8f} max {d:.3e}')
except Exception as e:
    print('baseline_regression failed:', type(e).__name__, e)