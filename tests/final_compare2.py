"""Final comparison: our corrected-parameter end-to-end output
(sample_data/dataset_HEP.set, saved WITH brainbeats by the re-attach patch)
vs the fixture reference (the exact MATLAB reference numbers).
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
import eegprep

ref = np.load(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python/tests/data/reference_hep_tf.npz',
              allow_pickle=True)
PY = eegprep.pop_loadset(
    r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
brb = PY.get('brainbeats')
if brb is None:
    print('brainbeats NOT saved in the file (re-attach failed)')
hp = (brb or {}).get('hrsp', {})
print('our hrsp keys:', sorted(hp.keys()) if isinstance(hp, dict) else hp)
hrsp_py = np.asarray(hp.get('hrsp'), float)
hepc_py = np.asarray(hp.get('hrpc', hp.get('hepc')), float)
hrsp_ref = np.asarray(ref['hrsp'], float)
hepc_ref = np.asarray(ref['hrpc'], float)
print('our hrsp shape:', hrsp_py.shape, ' ref:', hrsp_ref.shape)
if hrsp_py.shape == hrsp_ref.shape:
    r_b = np.corrcoef(hrsp_py.ravel(), hrsp_ref.ravel())[0, 1]
    r_c = np.corrcoef(hepc_py.ravel(), hepc_ref.ravel())[0, 1]
    d_b = np.abs(hrsp_py - hrsp_ref)
    print(f'END-TO-END HRSP corr {r_b:.4f}  max|dB| {d_b.max():.3f}  median {np.median(d_b):.4f}')
    print(f'END-TO-END HRPC corr {r_c:.4f}')
else:
    print('shape mismatch: different beat counts; corr over the common min shape')
    n = min(hrsp_py.shape[-1], hrsp_ref.shape[-1])
    m = min(hrsp_py.shape[0], hrsp_ref.shape[0])
    r_b = np.corrcoef(hrsp_py[:m, :, :n].ravel(), hrsp_ref[:m, :, :n].ravel())[0, 1]
    r_c = np.corrcoef(hepc_py[:m, :, :n].ravel(), hepc_ref[:m, :, :n].ravel())[0, 1]
    print(f'END-TO-END HRSP corr {r_b:.4f} (first {n} TF steps)'
          f'  HRPC corr {r_c:.4f}')
data = np.asarray(PY['data'], float)
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(PY['chanlocs']).ravel()]
roi = ['F1','F2','F3','F4','Fz','FC1','FC2','FC3','FC4','FCz',
       'C1','C2','C3','C4','Cz']
ridx = [i for i, l in enumerate(labels) if l in roi]
hep_py = data[ridx].mean(axis=0).mean(axis=1)
hep_ref = np.asarray(ref['hep_roi_from_epochs'], float)
n = min(hep_py.size, hep_ref.size)
r_h = np.corrcoef(hep_py[:n], hep_ref[:n])[0, 1]
print(f'ROI HEP epoch-mean corr {r_h:.4f} over first {n} pts')