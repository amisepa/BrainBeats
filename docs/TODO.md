# BrainBeats TODO

## Before the v1.6 release
- [ ] Run `brainbeats_tutorial.m` interactively once: plotting paths are not exercised by
      tests/run_tutorial_headless.m (IBI histogram, HEP plots, HRSP, features topos, coherence).
- [ ] Decide the release date in README "Version history" (written as 9/2026).
- [ ] Delete functions/corr_cmap.mat (unused).

## Next
- GUI options for 'hep_window', 'ppg_transit' and 'hep_baseline' (command line only for now);
  the GUI has no 'coherence' analysis; unused GUI options (RR-correction methods, legacy PPG
  detector options, eeg_interp) could be removed.
- HRSP is plotted for one channel only; export it for all channels (e.g. HEP.brainbeats.hrsp).
- Surrogate R-peaks (rigid shift, Park et al. 2018) as a HEP/HRSP null control.
- rm_heart removes ~40% of the heart-locked EEG amplitude on the sample data with one
  ICLabel heart component; a better CFA removal would need new methods (validate first).
- 'individualized' EEG bands: only alpha is individualized (per-channel peaks are unreliable
  for the other bands).
- get_RR fails on recordings shorter than ~3 s (filter length).
- Without an ECG, 'ppg_transit' needs the delay from the user: a default per sensor site
  (finger, ear, wrist) could be offered.
