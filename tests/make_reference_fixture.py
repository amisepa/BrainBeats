"""T1 parity fixture: HRSP/HRPC/HEP reference from dataset_HEP.set.

The MATLAB pipeline saved this file (sample_data/dataset_HEP.set from
brainbeats_process with analysis='hep'). Its brainbeats.hrsp struct carries the
reference HRSP (dB), HRPC ('hepc') and raw power, all 63 cleaned channels, plus
the epoched cleaned data (63 x 224 x 267) from which the ROI HEP reference is
recomputed (the HEP of the stored epochs == the HEP computed from the cleaned
continuous signal, because the epoch-level cleaning is linear).
"""
import os
import numpy as np
import eegprep

os.environ.setdefault('MPLBACKEND', 'Agg')
EEG = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
hrsp = EEG['brainbeats']['hrsp']
data = np.asarray(EEG['data'])
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(EEG['chanlocs']).ravel()]
roi = ['F1','F2','F3','F4','Fz','FC1','FC2','FC3','FC4','FCz',
       'C1','C2','C3','C4','Cz']
ridx = [i for i, l in enumerate(labels) if l in roi]
hep_ref = data[ridx].mean(axis=0).mean(axis=1)      # (224,)
out = {
    'hep_times': np.asarray(hrsp['hep_times']),
    'times': np.asarray(hrsp['times']),
    'freqs': np.asarray(hrsp['freqs']),
    'hrsp': np.asarray(hrsp['hrsp']),
    'hrpc': np.asarray(hrsp['hepc']),
    'power': np.asarray(hrsp['power']),
    'nBeats': int(hrsp['nBeats']),
    'hep_roi_from_epochs': hep_ref,
    'roi_channels': np.array(labels)[ridx],
    'all_channels': np.array(labels),
}
d = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                 'python', 'tests', 'data')
os.makedirs(d, exist_ok=True)
np.savez(os.path.join(d, 'reference_hep_tf.npz'), **out)
print('saved', {k: getattr(np.asarray(v), 'shape', v) for k, v in out.items()
                if k not in ('roi_channels', 'all_channels')},
      'nBeats', out['nBeats'], 'roi chans', len(ridx))