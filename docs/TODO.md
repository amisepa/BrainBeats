# BrainBeats TODO

## Before the v1.6 release
- [x] Update the copy of the plugin in eeglab/plugins (git pull) before testing in EEGLAB.
- [ ] Click through the new GUI: 'BrainBeats' in the EEGLAB menu bar, main window (logo,
      section headings, Load dataset / Sample dataset, analysis, heart channels and the '...'
      channel list, 'Perform analysis on'), parameters window (heart artifact removal: ICA,
      ECG regression or none). The window layouts and all plots have only been checked by
      automated tests that do not draw them.
- [x] Rerun `run_tutorial_headless(9)` ('rm_heart'): passes (160 s, 1 heart IC, 2.10 -> 1.12 uV,
      2026-10-07); the two timeouts were a busy machine.
- [x] Run `tests/test_gui_params` and the full `tests/run_tutorial_headless` once more: pass
      (26/26), with `tests/run_command_lines` (tutorial and help commands as written, 11/11;
      2026-10-07).
- [ ] Look at the 9 figures of the README (made by `tests/make_readme_figures`, which now runs
      headless too) in the MATLAB desktop once.
- [x] ASR cutoff default: 30 SD in the GUI and the command line (user, 2026-11-02).
- [ ] Decide the release date in README "Version history" (written as 9/2026).

## From the code review of v1.6 (fixed ones are in the README version history)
Scientific, may change results:
- Brain-heart coherence (experimental, out of the GUI): redesign at HRV time scales (NN series
  not high-passed at 0.5 Hz, MVAR at a matching rate); fit at the data rank with QR, not
  `(Y*Z')/(Z*Z')` and inv() in fdMVAR_5order, since cleaned EEG is rank-deficient.
- HRSP normalized over a window that includes the QRS: synthetic check (broadband transient at
  0 ms in noise, FDR points at 300-600 ms); if there is an offset, normalize by the mean power
  of the continuous data per channel and frequency.
- Surrogate beats can fall on epochs rejected by the cleaning (wider null, depends on the
  artifact load): pass the rejected spans to compute_hep_tf as exclusions.
- HRV normalized Lomb-Scargle: LF comes from 125-s windows and HF from 34-s windows, each
  normalized by its own variance, and band powers scale with the heart rate. Change the default
  to 'LombScargle' (ms^2), LF and HF from one window of 2 min or more, plus LF/(LF+HF)?
- clean_rr keeps ectopic pairs (0.75 x / 1.25 x the local interval pass both rules; RMSSD +40%
  at 30 ms over 4 min): injection check, then a successive-difference rule.
- PAT fallback (250 ms) is a warning: error unless 'ppg_transit' is given? Show the PAT method
  in the HEP figure titles; check the censoring of the pairing in 'peaks' mode.
- Validation: every number in the README comes from one subject. Run other ds003838 subjects,
  an annotated ECG/PPG benchmark and an HRV cross-check; v1.4 vs v1.6 table of changed outputs.
- Alpha asymmetry is ln(L) - ln(R), the opposite of the usual convention: flip it or state it.
- PRSA uses 20% (PhysioNet default); Bauer et al. (2006) use 5%.
Robustness:
- Absolute thresholds: rank `eig(cov) > 1e-7` (clean_eeg, run_HEP, apply_car),
  flatline `< 0.05` in get_RR (fails for ECG in volts), `signal_filt(fs:...)` with non-integer fs.
- PPG SQI clipping test (`|diff| > 0.5`) is in absolute units: normalize each beat.
- Fewer than 30 channels: no re-reference before ICLabel: allow average with a warning, or say so.
- Bad-epoch HF metric (data minus 45-Hz low-pass) does nothing after a 30-Hz low-pass, and the
  padded epochs overlap: evaluate on the analysis window with a band under the low-pass.
- IBIs are computed after synthetic beats are dropped: compute them before.
- 'rr' mode cannot take a PAT: accept a numeric 'ppg_transit'; tests/test_rr_vs_ppg.m locks HEPs
  to pulses.
- Nonlinear EEG features are O(N^2) on whole recordings and on by default: use 30-60 s windows.
- Generic helper names on the path: package them in the redesign.
- Tests: most headless cases only check that the run finishes: add tests/readme_numbers.m (the
  numbers quoted in the README, with tolerances) and a synthetic HRV and EEG features test.
Documentation:
- References to check: 'Park et al. (2018)' in compute_hep_tf.m and the README (surrogate
  heartbeats; Park et al., 2014, Nature Neuroscience?), Steinfath 2026, Ros 2025, Lee 2024; "used
  and cited in" mixes uses and reviews.
- Plugin header says GPL-2+, LICENSE is GPL-3; GEDAI is installed without asking.
- tests/prep_sampledata.m says 125 Hz; the sample dataset is 250 Hz.

## Redesign (after v1.6)
Goal: steps that can be combined, instead of 4 exclusive analyses; nothing slow computed
unless asked; heart-component removal as an optional preprocessing step (done for HEP in
v1.6: 'heart_removal'; the 'rm_heart' analysis was removed).

1. Steps, each a function usable alone and recorded in EEG.brainbeats:
   - heart: beat detection (ECG/PPG/rr), QA (SQI, rpeak_qa), clean_rr, PPG transit time;
     writes 'R-peak' events and EEG.brainbeats.heart (beats, NN, QA)
   - EEG cleaning: filter/reference/bad channels, ASR or bad epochs, ICA + ICLabel
   - heart-component removal (optional, before any analysis): reuse the cleaning ICA (no
     second ICA), flag heart ICs (ICLabel >= threshold), report the heart-locked residual
   - measures (any subset): HEP (window, baseline options), HRSP/HRPC, surrogate control,
     HRV features, EEG features, brain-heart coherence
2. brainbeats_process becomes a driver: 'steps' (or checkboxes in the GUI) select what
   runs; the current 'analysis' values map onto steps so existing scripts keep working
   ('hep' = heart + cleaning + HEP).
3. GUI in three parts: (1) data & heart signal, (2) preprocessing (EEG cleaning, heart-
   component removal), (3) measures: one checkbox per measure with an "Options" button,
   all off by default except the ones the user ticks. EEGLAB's inputgui has no tabs:
   either three successive inputgui windows (EEGLAB style, testable with the inputgui stub
   used for the v1.6 GUI checks) or one custom figure with uitabgroup tabs.
4. Outputs: the continuous (cleaned) dataset with EEG.brainbeats.{heart, preprocessing,
   features, coherence, parameters}, and the epoched HEP dataset with EEG.brainbeats.{hep,
   hrsp, surrogate} returned as a second output and saved as <filename>_HEP.set.
5. Tests: one headless test per step, plus the GUI stub checks for the new windows.

Decisions needed: tabs vs successive windows; second output for the HEP dataset; whether
heart-component removal keeps the heart channel in the ICA (as 'rm_heart' did) or
reuses the cleaning ICA (faster, one decomposition).

## Next (other)
- Heart artifact removal on the sample data: ICA with the ECG channel (former 'rm_heart', threshold
  .65) removes 47% of the heart-locked EEG amplitude and the ECG regression 72%, but the default
  HEP cleaning (ICA without the ECG channel, threshold .75) only 2% (2.18 -> 2.14 uV; 21% at 200-500
  ms). 6 ds003838 subjects (2026-11-02, tests/compare_heart_removal.m): 'ica' removed nothing on
  any; regression 29% of the QRS amplitude on average (3-70%) and kept 80% (54-100%) of a
  synthetic response; ICA with the ECG channel ('rm_heart', removed) took the response on 5/6.
  Default: 'ica' (user, 2026-11-02). Candidates from the literature (2026-11-02), to test with
  compare_heart_removal before adding any:
  - Heart ICs by correlation of the IC time course with the ECG (SASICA, threshold .6): sensitivity
    .975 vs .21 for ICLabel on 21,375 ICs (Champetier et al., CARACAS, bioRxiv 2025,
    doi:10.1101/2025.09.02.673728); same idea as MNE ICA.find_bads_ecg ('correlation' or 'ctps').
    Likely why ICLabel found no heart IC here. Risk: also removes heartbeat-locked brain activity.
  - CARACAS (same paper): heart ICs without an ECG (R-wave detection on each IC), in SASICA.
  - MNE compute_proj_ecg: SSP projectors from ECG-locked epochs (removes a fixed topography).
  - Controls rather than perfect removal: surrogate-heartbeat correction (Steinfath et al., 2026,
    Psychophysiology, doi:10.1111/psyp.70297), pseudotrial correction for task data (Imaging
    Neuroscience 2025, doi:10.1162/IMAG.a.30), surface Laplacian check, analysis windows after the
    T wave. Isoelectric-patient study (Kandasamy et al., bioRxiv 2026): ICA weakened the HEP, and
    the cardiac artefact topography is not predicted by the ECG.
  Open: ICA with the ECG channel only for components not locked to the heartbeat? A realistic
  synthetic HEP (variable amplitude and latency, cortical topography).
- Surrogates also for HEP contrasts between conditions (e.g. per-condition beat subsets).
- 'individualized' EEG bands: only alpha is individualized (per-channel peaks are unreliable
  for the other bands).
- get_RR fails on recordings shorter than ~3 s (filter length).
- PAT without an ECG (tests/compare_ppg_transit.m, 6 ds003838 subjects, 2026-11-02): the EEG
  estimate was wrong on 2/6 (628 and 116 ms vs 408 and 392 ms from the ECG) with ratios 2.2/2.1,
  like a correct one (2.1); the 250 ms default was wrong on 6/6 (true 384-424 ms). With an ECG the
  correction works (r .52-.86 with the ECG-locked HEP). Decide: no silent correction without ECG?
- PPG polarity: the ds003838 PPG (sample included) is in light-intensity convention (fast DOWNWARD
  edge at each pulse), so 'valleys' marks the systolic volume peak (R + 384-424 ms) and 'peaks' a
  point just before the pulse onset (R + 108-128 ms, smoothed). 'peaks' matched the ECG-locked HEP
  better on 4/6 subjects (beat-to-beat IQR 9-24 vs 16-44 ms). Detect the polarity (as for the ECG)
  so that 'valleys' is the pulse onset whatever the device; then locate the onset on a less
  smoothed signal.
- 'keep_heart' stores the raw heart channel (sample ECG has a large DC offset): filter it?
