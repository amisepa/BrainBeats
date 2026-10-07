"""Fix roi labels + .beat extraction on the regenerated reference."""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import scipy.io as sio

p = r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set'
with open(p, 'rb') as fh:
    hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
bb = hdr['brainbeats']
roi = bb.roi
rf = getattr(roi, '_fieldnames', None)
print('roi type:', type(roi), 'fields:', rf)
if rf:
    print('roi.channels?', [f for f in rf])
    for f in rf:
        v = getattr(roi, f)
        try:
            print(' ', f, np.asarray(v, str).ravel()[:4])
        except Exception:
            print(' ', f, type(v))
else:
    roi = bb.roi = bb.roi
# epoch events
ep = hdr['epoch']
print('epoch type:', type(ep), 'size:', ep.size if hasattr(ep, 'size') else None)
e0 = ep[0] if isinstance(ep, np.ndarray) else ep
ev = e0.event
print('event type:', type(ev))
ef = getattr(ev, '_fieldnames', None)
print('event fields:', ef)
if ef:
    for f in ef:
        v = getattr(ev, f)
        try:
            print(' ', f, np.ravel(v)[:3])
        except Exception:
            print(' ', f, type(v))
# maybe events live in EEG.event (single events list) with .beat
EEGev = getattr(hdr.get('EEG'), 'event', None) if hasattr(hdr.get('EEG', None), '_fieldnames') else None
print('EEG wrapper present:', type(hdr.get('EEG')))
# top-level event key?
top_ev = hdr.get('event')
print('top-level event:', type(top_ev),
      (getattr(top_ev, '_fieldnames', None) if hasattr(top_ev, '_fieldnames')
       else (np.asarray(top_ev, dtype=object).shape if top_ev is not None else None)))