"""Final regen fixture: roi struct (channels/hep/hep_times/tf) + .beat from events.

The stored events are a structured ARRAY (fields via dtype names), and roi is a
MAT struct with channels + hep + hep_times + tf -- the full ROI reference.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import scipy.io as sio

p = r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set'
with open(p, 'rb') as fh:
    hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
bb = hdr['brainbeats']
roi = bb.roi
out = dict(np.load(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python/tests/data/regen_reference.npz',
                   allow_pickle=True))
out['roi_channels'] = np.asarray(getattr(roi, 'channels'), str).ravel()
out['roi_hep'] = np.asarray(getattr(roi, 'hep'), float)
out['roi_hep_times'] = np.asarray(getattr(roi, 'hep_times'), float).ravel()
top_ev = hdr['event']
names = top_ev.dtype.names
print('event fields:', names)
if names and 'beat' in names:
    beats = np.asarray([np.ravel(e['beat'])[0] for e in top_ev], float)
    out['epoch_beats'] = beats
    print('beats n=', beats.size, 'range', beats.min(), beats.max())
dst = r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python/tests/data/regen_reference.npz'
np.savez(dst, **out)
print('saved', dst)
print('roi_channels:', out['roi_channels'])
print('roi_hep:', out['roi_hep'].shape)