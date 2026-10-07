"""Standalone verification of the brainbeats re-attach on the last saved .set:
does loadmat see EEG as an object array? Does injecting brainbeats survive?
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import scipy.io as sio

p = r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set'
with open(p, 'rb') as fh:
    hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
print('top keys:', sorted(hdr.keys()))
e = hdr.get('EEG')
print('EEG type:', type(e))
if isinstance(e, np.ndarray):
    print('  dtype:', e.dtype, 'shape', e.shape)
    inner = e[0]
    print('  inner dtype', type(inner), 'fields' if not isinstance(inner, dict) else dict(keys))
    if hasattr(inner, '_fieldnames'):
        print('  fields:', inner._fieldnames)
elif e is not None:
    print('  fields:', getattr(e, '_fieldnames', None))
# now test a round-trip: inject a dict into EEG.brainbeats and save
test = {'hrsp': {'probe': np.arange(6.0)}}
if hasattr(e, '__dict__'):
    e.brainbeats = test
    sio.savemat(p + '.probe.set', hdr, oned_as='row')
    back = sio.loadmat(p + '.probe.set', struct_as_record=False, squeeze_me=True)
    print('roundtrip brainbeats:', type(back['EEG'].brainbeats),
          back['EEG'].brainbeats if isinstance(back['EEG'].brainbeats, dict)
          else getattr(back['EEG'].brainbeats, 'hrsp', 'MISSING'))
else:
    print('EEG is an object array: inject via item assignment before saving')
os.remove(p + '.probe.set') if os.path.exists(p + '.probe.set') else None