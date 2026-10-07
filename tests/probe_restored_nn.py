"""Does the git-restored reference carry the beat train (NN/NN_times)?"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import scipy.io as sio

p = r'C:/Users/ccann/DATASET_HEP_reference.set'
with open(p, 'rb') as fh:
    hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
bb = hdr.get('brainbeats')
if bb is None:
    print('no brainbeats')
else:
    pre = getattr(bb, 'preprocessings', None)
    fields = getattr(pre, '_fieldnames', None)
    print('preprocessings fields:', fields)
    if fields:
        nnt = getattr(pre, 'NN_times', None)
        if nnt is not None:
            nnt = np.asarray(nnt, float).ravel()
            print('NN_times n =', nnt.size, 'first', nnt[:4])
            np.save(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/tests/'
                    r'nn_times_restored.npy', nnt)
            print('SAVED nn_times_restored.npy')
    else:
        print('preprocessings:', pre)