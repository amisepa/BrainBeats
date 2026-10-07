"""Verify the restored MATLAB reference .set and finalize the parity result.

Compares our end-to-end output (still in the overwritten sample_data file!)
against the restored reference. The overwritten file at
sample_data/dataset_HEP.set is OUR run_hep output (Python). The restored
reference lives at C:/Users/ccann/DATASET_HEP_reference.set.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
import eegprep

REF = eegprep.pop_loadset(r'C:/Users/ccann/DATASET_HEP_reference.set')
PY = eegprep.pop_loadset(
    r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
print('reference:', np.asarray(REF['data']).shape,
      'python-run:', np.asarray(PY['data']).shape)
br = REF.get('brainbeats', {})
print('reference brainbeats keys:', sorted(br.keys()) if isinstance(br, dict) else 'none')
hrsp_ref = br.get('hrsp', {})
if isinstance(hrsp_ref, dict) and 'hrsp' in hrsp_ref:
    b = PY['brainbeats']['hrsp']
    hrsp_ref = np.asarray(br['hrsp']['hrsp'], float)
    hepc_ref = np.asarray(br['hrsp']['hepc'], float)
    hp = PY['brainbeats'].get('hrsp', {})
    hrsp_py = np.asarray(hp.get('hrsp'), float)
    hepc_py = np.asarray(hp.get('hrpc', hp.get('hepc')), float)
    if hrsp_py.shape == hrsp_ref.shape:
        r_b = np.corrcoef(hrsp_py.ravel(), hrsp_ref.ravel())[0, 1]
        r_c = np.corrcoef(hepc_py.ravel(), hepc_ref.ravel())[0, 1]
        print(f'END-TO-END (Picard) HRSP corr {r_b:.4f}  HRPC corr {r_c:.4f}')
    else:
        print('shapes differ:', hrsp_py.shape, hrsp_ref.shape)
print('reference kept trials:', REF['trials'], ' our trials:', PY['trials'])