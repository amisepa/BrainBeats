# BrainBeats tests and validation scripts

Each script starts with a help block (purpose, usage, and for the validation scripts the
results and the date they were obtained). Run them from this folder, with EEGLAB on the path
(`eeglab; close`).

## Before a release
| Script | What it checks | Time |
|---|---|---|
| `run_tutorial_headless.m` | Every analysis of the tutorial and the options added since (32 cases), plots off; `run_tutorial_headless(ids)` runs a subset | 30-60 min |
| `run_command_lines.m` | The tutorial and the help examples, run as they are written | ~30 min |
| `test_gui_params.m` | GUI defaults and mappings for every analysis without opening windows, and one GUI-driven run with its history command | ~5 min |
| `test_hep_figures.m` | Content of the HEP figures (drawn invisible) | ~5 min |
| `test_eeglab_menu.m` | The 'BrainBeats' menu of EEGLAB stays usable after a HEP run (opens EEGLAB) | ~5 min |
| `test_compute_hep_tf.m` | HEP/HRSP/HRPC and surrogate statistics on synthetic data (no EEGLAB needed) | ~2 min |
| `test_clean_rr_ectopic.m` | Ectopic-beat replacement and RR_ALTERNATION on synthetic RR series | seconds |
| `test_thresholds.m` | Thresholds that must not depend on the data's units (rank, ECG flatline, PPG clipping) | seconds |
| `test_fuzzy_fractal.m` | Fuzzy entropy and fractal dimension from ASCENT (values, warnings, NaN without ASCENT) | seconds |
| `test_features_synthetic.m` | HRV and EEG features against synthetic data with known values, and through the whole 'features' pipeline from a synthetic ECG | ~1 min |
| `make_readme_figures.m` | Regenerates the README figures (`figures/`) from the sample data | ~15 min |

## Validation (numbers quoted in the README and in function help)
Most need the ds003838 recordings prepared by `prep_ds003838.m` (OpenNeuro, not in the repo).

| Script | Question |
|---|---|
| `validate_ectopic_mitbih.m` | Ectopic rules on the MIT-BIH Arrhythmia Database (download from PhysioNet) |
| `validate_ppg_polarity.m` | PPG polarity statistics on 70 PPGs of known polarity, clean and with injected artifacts |
| `validate_subsample_rr.m` | Sub-sample beat times vs the sample grid: synthetic ECG, ds003838 ECG and PPG, MIT-BIH, resampled to 128-250 Hz |
| `check_rr_alternation.m` | RR_ALTERNATION on healthy resting recordings (no false warning) |
| `compare_ectopic_rules.m` | Ectopic rules on real RR series with injected premature, missed and extra beats |
| `compare_hrv_windows.m` | HRV band powers with LF and HF from the same 120-s windows |
| `compare_badepoch_metrics.m` | Bad-epoch detection for HEP, with injected muscle bursts |
| `compare_heart_removal.m` | Heart artifact left and synthetic brain response kept by each removal method |
| `compare_regression_copy.m` | ECG regression applied to the ICA copy too, or not |
| `compare_iclabel_source.m` | ICLabel on the broadband ICA copy vs on the 0.5-30 Hz data |
| `compare_ppg_transit.m` | PPG pulse arrival time correction vs the ECG-locked HEP |
| `ica_compare.m` | ICA settings of the HEP cleaning |
| `sim_hep_fwer.m` | False positives and power of the HEP surrogate statistics |
| `sim_surrogate_reference.m`, `check_surrogate_reference.m` | Reference of the surrogate test (simulation, real recordings) |
| `sim_cfa_drift.m` | Sensitivity of the QRS-map change test |

## Beats detected elsewhere ('rr' mode), Muse 2 data in `data/`
`check_alignment.m`, `check_drift.m` (timestamps of the Muse EEG, optics and beat
detections), `test_clean_rr.m` (clean_rr on beats already detected), `test_rr_vs_ppg.m` (HEP
from the Muse beats vs BrainBeats PPG detection; figures in `figures/`).

## Helpers
`highpass_eeg.m` (high-passed copy of a test dataset), `prep_ds003838.m` (the ds003838
recordings used above), `prep_sampledata.m` (how `sample_data/dataset.set` was built).
