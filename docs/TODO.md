# BrainBeats TODO

## Before the v1.6 release
- [ ] Run `brainbeats_tutorial.m` interactively once: plotting paths are not exercised by
      tests/run_tutorial_headless.m (IBI histogram, HEP plots, HRSP/HEPC, surrogate HEP,
      features topos, coherence).
- [ ] Decide the release date in README "Version history" (written as 9/2026).

## Redesign (after v1.6)
Goal: steps that can be combined, instead of 4 exclusive analyses; nothing slow computed
unless asked; heart-component removal as an optional preprocessing step.

1. Steps, each a function usable alone and recorded in EEG.brainbeats:
   - heart: beat detection (ECG/PPG/rr), QA (SQI, rpeak_qa), clean_rr, PPG transit time;
     writes 'R-peak' events and EEG.brainbeats.heart (beats, NN, QA)
   - EEG cleaning: filter/reference/bad channels, ASR or bad epochs, ICA + ICLabel
   - heart-component removal (optional, before any analysis): reuse the cleaning ICA (no
     second ICA), flag heart ICs (ICLabel >= threshold), report the heart-locked residual
   - measures (any subset): HEP (window, baseline options), HRSP/HEPC, surrogate control,
     HRV features, EEG features, brain-heart coherence
2. brainbeats_process becomes a driver: 'steps' (or checkboxes in the GUI) select what
   runs; the current 'analysis' values map onto steps so existing scripts keep working
   ('hep' = heart + cleaning + HEP; 'rm_heart' = heart + cleaning + heart removal).
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
heart-component removal keeps the heart channel in the ICA (as rm_heart does now) or
reuses the cleaning ICA (faster, one decomposition).

## Next (other)
- The GUI has no 'coherence' analysis; unused GUI options (RR-correction methods, legacy PPG
  detector options, eeg_interp) could be removed (part of the redesign).
- rm_heart removes ~40-50% of the heart-locked EEG amplitude on the sample data with one
  ICLabel heart component; a better cardiac field removal needs new methods (validate first).
- Surrogates also for HEP contrasts between conditions (e.g. per-condition beat subsets).
- 'individualized' EEG bands: only alpha is individualized (per-channel peaks are unreliable
  for the other bands).
- get_RR fails on recordings shorter than ~3 s (filter length).
- Without an ECG, 'ppg_transit' needs the delay from the user: a default per sensor site
  (finger, ear, wrist) could be offered.
