"""
brainbeats-py -- Python port of the BrainBeats EEGLAB plugin, HEP/HRSP/HRPC mode.

This package is an eegprep extension: it adds the BrainBeats heart analyses
(heartbeat-evoked potentials, heartbeat-related spectral perturbation and
heartbeat-evoked phase coupling) to the Python EEGLAB (eegprep, sccn).

Layout mirrors the MATLAB plugin:

    functions/          one module per MATLAB function (same names, snake_case)
    tests/              parity and unit tests against the MATLAB reference output

Validated building blocks reused from the HEP_neurofeedback project (do NOT
reimplement here -- they are imported from their canonical file):
    get_rr, resolve_polarity   <- get_rr_v25.py     (validated vs 277,320 reference peaks)
    clean_rr, qrs_bandpass     <- rr_cleaning_v25.py (validated vs 70 reference files, 100.0000%)
"""
from __future__ import annotations

__version__ = '0.1.0'
from .csd_transform import csd_transform_data, csd_matrix, get_gh  # noqa: E401
from .estimate_pat import estimate_pat  # noqa: E401
from .estimate_pat_eeg import estimate_pat_eeg  # noqa: E401
from .remove_heart_regression import remove_heart_regression, cfa_amplitude  # noqa: E401
from .compute_fe import compute_fe  # noqa: E401
from .fractal_volatility import fractal_volatility  # noqa: E401
from .resample_NN import resample_NN  # noqa: E401
from .compute_psd import compute_psd, nextpow2  # noqa: E401
from .resting_iaf import resting_iaf, get_freqBounds, findF1, findF2  # noqa: E401
from .compute_asymmetry import compute_asymmetry  # noqa: E401
from .get_hrv_features import get_hrv_features  # noqa: E401
from .get_eeg_features import get_eeg_features  # noqa: E401
from .extract_features import extract_features  # noqa: E401
from .run_features import run_features, pop_run_features  # noqa: E401
