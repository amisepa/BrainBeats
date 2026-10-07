"""Identify the current dataset_HEP.set: whose output is it, and does it carry
brainbeats (MATLAB form or our re-attached form)?
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import scipy.io as sio

p = r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set'
with open(p, 'rb') as fh:
    hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
bb = hdr.get('brainbeats')
pre = getattr(bb, 'preprocessings', None) if bb is not None else None
if pre is not None and hasattr(pre, '_fieldnames'):
    nnt = getattr(pre, 'NN_times', None)
    comp = getattr(pre, 'removed_eeg_components', None)
    print('PREPROCESSINGS PRESENT (MATLAB structure or our re-attach)')
    print(' fields:', pre._fieldnames[:8])
    if nnt is not None:
        print(' NN_times n:', np.asarray(nnt).size)
    if comp is not None:
        print(' removed comps:', np.ravel(comp))
else:
    print('NO brainbeats preprocessings -> ', 'brainbeats None' if bb is None else type(bb))
print('nbchan/pnts/trials:', hdr.get('nbchan'), hdr.get('pnts'), hdr.get('trials'))
w = hdr.get('icaweights')
print('icaweights shape:', np.asarray(w).shape if w is not None else None)
print('setname:', hdr.get('setname'))