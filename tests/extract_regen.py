"""Extract the regenerated MATLAB PPG reference into a parity fixture.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import scipy.io as sio

p = r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set'
with open(p, 'rb') as fh:
    hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
bb = hdr['brainbeats']
pre = bb.preprocessings
fields = getattr(pre, '_fieldnames', [])
print('preprocessings fields:', fields)
nnt = np.asarray(getattr(pre, 'NN_times'), float).ravel()
nn = np.asarray(getattr(pre, 'NN'), float).ravel()
hrsp = bb.hrsp
hf = getattr(hrsp, '_fieldnames', [])
print('hrsp fields:', hf)
out = {
    'NN_times': nnt, 'NN': nn,
    'hrsp': np.asarray(getattr(hrsp, 'hrsp'), float),
    'hrpc': np.asarray(getattr(hrsp, 'hrpc'), float),
    'power': np.asarray(getattr(hrsp, 'power'), float),
    'times': np.asarray(getattr(hrsp, 'times'), float).ravel(),
    'freqs': np.asarray(getattr(hrsp, 'freqs'), float).ravel(),
    'nBeats': np.asarray(getattr(hrsp, 'nBeats'), float).ravel(),
    'channels': np.asarray(getattr(hrsp, 'channels'), str).ravel(),
    'data_shape': np.asarray(hdr['EEG'].data if hasattr(hdr.get('EEG', None), 'data')
                             else hdr['data'], float).shape,
}
# stored epoch .beat values (event field of each epoch)
ep = hdr['epoch']
beats = []
for e in range(ep.size):
    ev = ep[e].event
    b = getattr(ev, 'beat', None) if hasattr(ev, '_fieldnames') else None
    if b is not None:
        beats.append(float(np.ravel(b)[0]))
print('epochs with beat:', len(beats), '/', ep.size)
out['epoch_beats'] = np.asarray(beats, float)
roi = getattr(bb, 'roi', None)
if roi is not None:
    out['roi'] = np.asarray(roi, str).ravel()
dst = r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python/tests/data/regen_reference.npz'
np.savez(dst, **out)
print('saved', dst)
print('nBeats:', out['nBeats'], ' NN:', nnt.size, ' hrsp shape:', out['hrsp'].shape)
print('roi:', out.get('roi'))