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