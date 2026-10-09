# BrainBeats TODO

## Before the v1.6 release
- [ ] Update the copy of the plugin in eeglab/plugins (git pull) before testing in EEGLAB.
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
- [ ] Run `tests/make_readme_figures` in the MATLAB desktop (plotting paths are not
      exercised headless): it saves the 11 figures README.md points to (gui_main,
      heartbeats_ecg, hep_ecg, hep_tf, hep_ppg, features_psd, features_eeg,
      rm_heart_components, rm_heart_data, coherence_spectra, coherence_bands). Check each,
      then `git rm` the unreferenced old figures (fig4, fig11, fig17, fig21, fig22, fig27,
      coherence_allfreqs, coherence_topo, diagram, logo2-4).
- [ ] GUI ASR cutoff default is 30, command line default 50: pick one.
- [ ] Decide the release date in README "Version history" (written as 9/2026).

## Redesign (after v1.6)
Goal: steps that can be combined, instead of 4 exclusive analyses; nothing slow computed
unless asked; heart-component removal as an optional preprocessing step (done for HEP in
v1.6: 'heart_removal'; the 'rm_heart' analysis remains in the command line).

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
- (done in v1.6: coherence in the GUI, unused GUI options removed; the command-line parser
  still accepts the legacy RR/PPG options, unused.) Old note: unused GUI options (RR-correction methods, legacy PPG
  detector options, eeg_interp) could be removed (part of the redesign).
- Heart artifact removal on the sample data: ICA removes ~40-50% of the heart-locked EEG
  amplitude (one ICLabel heart component), the ECG regression 72%. The regression was checked
  on one subject and one synthetic response (91% kept): validate it on more data.
- Surrogates also for HEP contrasts between conditions (e.g. per-condition beat subsets).
- 'individualized' EEG bands: only alpha is individualized (per-channel peaks are unreliable
  for the other bands).
- get_RR fails on recordings shorter than ~3 s (filter length).
- PAT from the EEG cardiac field artifact (ESTIMATE_PAT_EEG) is validated on one subject only
  (420 vs 424 ms with the ECG): check it on datasets with ECG + PPG + EEG (e.g. Heartbeam), and
  whether the 2x peak/median GFP threshold rejects bad estimates.
- 'keep_heart' stores the raw heart channel (sample ECG has a large DC offset): filter it?

## Python port (branch python-port)
- 'ppg_transit' values are in ms, not samples: ESTIMATE_PAT pairs each PPG beat with the
  closest preceding R-peak within 50-600 ms (no cross-correlation), and BRAINBEATS_PROCESS
  shifts the beats by round(pat/1000*srate). The 424 on the sample data is 424 ms (106
  samples at 250 Hz; 420 ms from the EEG alone), a normal pulse arrival time for a pulse
  onset (literature: ~250 ms at the finger, 200-450 ms overall), not 1.7 s. The port should
  read it in ms.
- The MATLAB reference changed in v1.6 after bd71e56 (PPG detection, HRPC name, 0.5-Hz HEP
  high-pass, surrogates, 'hep_level', 'heart_removal'): regenerate the reference fixtures.

