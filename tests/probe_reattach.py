"""Why does pop_loadset fail on the re-attached dataset_HEP.set? Minimal probe.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import scipy.io as sio

p = r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set'
try:
    with open(p, 'rb') as fh:
        hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
    print('top keys:', sorted(hdr.keys())[:8], '... total', len(hdr))
    print('brainbeats type:', type(hdr.get('brainbeats')))
    bb = hdr['brainbeats']
    if hasattr(bb, '_fieldnames'):
        print('brainbeats fields:', bb._fieldnames)
    elif isinstance(bb, dict):
        print('brainbeats keys:', sorted(bb.keys()))
except Exception as e:
    print('loadmat FAILED:', type(e).__name__, e)
# pop_loadset path
import eegprep
try:
    EEG = eegprep.pop_loadset(p)
    print('pop_loadset OK', np.asarray(EEG['data']).shape)
except Exception as e:
    print('pop_loadset FAILED:', type(e).__name__, str(e)[:200])