"""Verify the re-attached brainbeats roundtrips through BOTH loaders."""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import scipy.io as sio
import eegprep

p = r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set'
with open(p, 'rb') as fh:
    hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
bb = hdr.get('brainbeats')
print('brainbeats present:', bb is not None)
bf = getattr(bb, '_fieldnames', None) or (sorted(bb.keys()) if isinstance(bb, dict) else type(bb))
print('fields:', bf)
hrsp = getattr(bb, 'hrsp', None) if bf else None
if isinstance(bb, dict):
    hrsp = bb.get('hrsp')
if hrsp is not None:
    print('hrsp fields:', getattr(hrsp, '_fieldnames', None) or sorted(hrsp.keys())[:8])
EEG = eegprep.pop_loadset(p)
print('pop_loadset: brainbeats present:', 'brainbeats' in EEG)
print('pop_loadset: data', np.asarray(EEG['data']).shape)