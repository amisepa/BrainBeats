"""Roundtrip test on a known-good pop_saveset output: loadmat -> savemat -> loadmat.
Find the breaking field.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import scipy.io as sio

SRC = (r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/'
       r'dataset_HEP.set')          # CORRUPTED (post-reattach)
GOOD = (r'C:/Users/ccann/DATASET_HEP_reference.set')  # git-restored, valid MATLAB file

with open(GOOD, 'rb') as fh:
    hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
out = r'C:/Users/ccann/AppData/Local/hermes/profiles/superscientist/cache/scratch/rt_test.set'
sio.savemat(out, hdr, oned_as='row')
try:
    with open(out, 'rb') as fh:
        h2 = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
    print('roundtrip OK; keys', len(h2))
except Exception as ex:
    print('roundtrip FAILED:', type(ex).__name__, str(ex)[:160])
# field-by-field: which type breaks it
ee = hdr.get('EEG')
print('reference file top keys:', sorted(hdr.keys())[:8])
print('EEG wrapper?', type(ee))
os.remove(out)
os.remove(SRC + '.probe.set') if os.path.exists(SRC + '.probe.set') else None