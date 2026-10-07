"""Save our cleaned 1-based beat train for the stage-B alignment check."""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/HEP_neurofeedback/python')
import eegprep
from get_rr_v25 import get_rr
from rr_cleaning_v25 import clean_rr, qrs_bandpass

EEG = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset.set')
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(EEG['chanlocs']).ravel()]
ecg = np.asarray(EEG['data'])[labels.index('ECG')].astype(float)
fs = EEG['srate']
peaks_full, info = get_rr(ecg, fs, drop_first=False)
filt = qrs_bandpass(ecg, fs)
peaks_oneb = peaks_full[1:]
rr = np.diff(peaks_full) / fs
nn, npk, idx_bad, info_rr = clean_rr(
    rr, peaks_oneb, fs, peak_amp=filt[peaks_oneb - 1], ecg_signal=ecg,
    sig_t=np.arange(ecg.size) / fs, interpolate_missing=True)
npk = np.asarray(npk, float)
real = npk[np.isfinite(npk)]
np.save(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/tests/our_train.npy',
        real.astype(np.int64))
# HEP set events: beat numbering + epoch-relative latencies
HEP = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
beats = np.array([ev.get('beat') for ev in HEP['event']], float)
lat = np.array([ev.get('latency') for ev in HEP['event']], float)
print('HEP events: beat range', beats.min(), beats.max(),
      'consecutive:', bool(np.all(np.diff(beats) == 1)))
print('n events', len(beats))
print('our train n =', real.size, 'first 5:', real[:5], 'last 3:', real[-3:])
# Map: reference kept beat k (k = 1..301) has sample index = our kept train position IF the
# trains align. Check consistency: kept beats are a SUBSET of the clean train at the same
# positions (both from get_RR + clean_rr). Take reference beat indices (1-based into the
# MATLAB kept train) and test whether our train's beats at those same positions reproduce
# plausible latencies (increasing, within 227 s).
print('first HEP event beat number:', beats.min(), '-> reference removed beats 1-4 at the start.')
print('our train[4] (beat 5):', real[4], 'vs the sample artifact period ~0-3 s (750 samples):')