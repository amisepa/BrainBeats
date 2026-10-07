"""Extract MATLAB's ICA solution + component removal list from the regen reference.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import scipy.io as sio

p = r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set'
with open(p, 'rb') as fh:
    hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
w = np.asarray(hdr['icaweights'], float)
inv = np.asarray(hdr['icawinv'], float)
print('icaweights', w.shape, 'icawinv', inv.shape)
pre = hdr['brainbeats'].preprocessings
comp = np.asarray(getattr(pre, 'removed_eeg_components'), float).ravel()
print('removed_eeg_components:', comp)
np.savez(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python/tests/data/'
         r'regen_ica.npz', icaweights=w, icawinv=inv, removed=comp)
print('saved regen_ica.npz')