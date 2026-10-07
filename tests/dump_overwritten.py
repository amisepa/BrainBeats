"""Diagnose what the overwritten file actually contains now."""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import eegprep
PY = eegprep.pop_loadset(
    r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
print('keys:', sorted(PY.keys()))
for k in PY.keys():
    v = PY[k]
    if isinstance(v, dict):
        print(' ', k, '->', sorted(v.keys())[:10])
    elif isinstance(v, np.ndarray):
        print(' ', k, '-> array', v.shape)
    else:
        print(' ', k, '->', type(v).__name__, v if isinstance(v, (int, float, str)) else '')