"""Which stored struct holds the reference numbers, and what does the stored
epoch mean equal? The run_HEP.m chain is:
  Xc (continuous, no baseline reg) -> brainbeats.hrsp (all-ch hrsp/hrpc)
                                    -> brainbeats.roi.hep (ROI HEP)
  HEP.data (epochs, WITH baseline_regression) -> saved .set
So the stored EPOCH mean is the REGRESSED HEP; brainbeats.roi.hep is the
UNREGRESSED one. Verify both, then update the parity fixture accordingly.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
import eegprep

HEP = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
br = HEP['brainbeats']
data = np.asarray(HEP['data'], float)
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(HEP['chanlocs']).ravel()]
roi = ['F1','F2','F3','F4','Fz','FC1','FC2','FC3','FC4','FCz',
       'C1','C2','C3','C4','Cz']
ridx = [i for i, l in enumerate(labels) if l in roi]
hep_hand = data[ridx].mean(axis=0).mean(axis=1)      # stored-epoch mean (regressed)

print('brainbeats keys:', sorted(br.keys()))
roi_st = br.get('roi', None)
if isinstance(roi_st, dict):
    print('roi keys:', sorted(roi_st.keys()))
    hp = np.asarray(roi_st['hep'], float)
    if hp.ndim > 1:
        hp = hp[0]
    print('roi.hep n =', hp.size)
    print(f'corr(stored-epoch-mean, roi.hep): '
          f'{np.corrcoef(hep_hand, hp[:hep_hand.size])[0, 1]:.6f}')
hr = br['hrsp']
print('hrsp keys:', sorted(hr.keys()))
# any 'hep' in hrsp? (surrogate branch stores surr.hep.real from tf.hep)
for kk in ('hep', 'hep_times', 'channels', 'nBeats', 'times', 'freqs'):
    if kk in hr:
        v = hr[kk]
        try:
            print('hrsp.' + kk, np.asarray(v, float).shape if not isinstance(v, (str, list)) else type(v))
        except Exception:
            print('hrsp.' + kk, type(v))
hepc = np.asarray(hr['hepc'], float) if 'hepc' in hr else None
print('hrsp.hepc' if hepc is not None else 'no hepc',
      hepc.shape if hepc is not None else '')