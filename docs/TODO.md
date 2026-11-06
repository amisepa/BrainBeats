# BrainBeats TODO

## Before the v1.6 release
- [x] Update the copy of the plugin in eeglab/plugins (git pull) before testing in EEGLAB.
- [ ] Click through the new GUI: 'BrainBeats' in the EEGLAB menu bar, main window (logo,
      section headings, Load dataset / Sample dataset, analysis, heart channels and the '...'
      channel list, 'Perform analysis on'), parameters window (heart artifact removal: none,
      ICA, CARACAS, and with an ECG: ECG correlation, ECG regression, SSP; PPG pulse arrival
      time options). The window layouts and all plots have only been checked by automated
      tests that do not draw them.
- [x] Rerun `run_tutorial_headless(9)` ('rm_heart'): passes (160 s, 1 heart IC, 2.10 -> 1.12 uV,
      2026-10-07); the two timeouts were a busy machine.
- [x] Run `tests/test_gui_params` and the full `tests/run_tutorial_headless` once more: pass
      (26/26), with `tests/run_command_lines` (tutorial and help commands as written, 11/11;
      2026-10-07).
- [ ] Look at the 9 figures of the README (made by `tests/make_readme_figures`, which now runs
      headless too) in the MATLAB desktop once.
- [x] ASR cutoff default: 30 SD in the GUI and the command line (user, 2026-11-02).
- [x] Release date in README "Version history": 11/2026 (user asked to push v1.6, 2026-11-03).

## Next (user, 2026-11-03)
- [x] ICA quality for HEP (user: improve the decomposition rather than lower the ICLabel eye
      threshold, .90 kept). Done 2026-11-03 (6 recordings x 10 settings): new HEP default = continuous 1-Hz copy, no low-pass (line noise
      notched), ASR 80 ('ica_asr'), ICLabel on the copy. Picard vs Infomax: no consistent
      difference (method default left at Infomax). Open: Picard diverged on 2 of 30 runs (fell
      back to Infomax, 1300-3800 s). ICLabel on the broadband copy vs on the 0.5-30 Hz data
      (same weights; tests/compare_iclabel_source.m, 6 recordings, 2026-11-05): same ICs removed
      on 5, one more eye IC removed on the copy on 1; more ICs muscle-labeled on the copy but
      none above the .99 threshold. Kept on the copy.
- [x] Surrogate statistics more conservative than per-latency FDR (user): TFCE (default) and
      max statistic, family-wise, 1000 surrogate HEPs ('hep_surr_stat'); simulation in
      tests/sim_hep_fwer.m, numbers in the README.
- [ ] IC figure: fewer plots. Clustering within one recording would merge independent sources;
      rather sort by heartbeat locking (max |z|) and optionally show only ICs that pass the
      surrogate test or are classified brain. Across subjects, cluster ICs in an EEGLAB STUDY
      (scalp maps + HEP).
- [x] README: hep_ecg and hep_ppg replaced by one overlay (hep_ecg_ppg.png: ROI HEP from the
      ECG, the corrected and the uncorrected PPG).
- [ ] plot_NN labels the ECG in uV; ds003838 stores the ECG in raw amplifier units (50,000-
      730,000 at the 99.5th percentile across 21 recordings, the tutorial dataset included).
      Label from EEG.chanlocs/units when known, else 'a.u.'.
- [x] More ds003838 subjects: all 65 usable screened (2026-11-03; 064, 065 fail with < 10
      heartbeats, not investigated). sub-047 stays the sample; sub-062 is the alternative.
- [x] HEP cleaning, ASR on the ICA copy (user, 2026-11-03: lax, 80 SD): default (item above).
      ASR 80 changed 0-0.4% of the samples of the ds003838 copies, 1.5-2.9% of the tutorial
      dataset's (added artifacts). ICA's rank: relative tolerance since 2026-11-05 (data_rank).

## Next (user, 2026-11-04)
- [x] Pre-R-peak values of each heartbeat exported for group covariates and single-trial
      models (HEP.brainbeats.prestim, <file>_HEP_beats.csv): baseline voltage (Alday's
      regressor) and induced theta/alpha/beta power. No correction inside BrainBeats.
- [x] Recording-level Alday correction removed (user, 2026-11-04): 'hep_baseline' only warns
      (BrainBeats:removed); 'hep_baseline_win' sets the exported baseline window;
      BASELINE_REGRESSION kept as a function for epochs that pool conditions.
- [x] QRS map stored ('hep_cfa', default on; CFA_MAP) with notes; PNG and .txt with 'save';
      README note on posture and head orientation (Dirlich et al. 1997; Plantard et al. 2025).
- [x] Within-recording change warning: sees only large lasting changes (tests/sim_cfa_drift.m:
      no false alarm on 24 recordings; 30-deg turn 3-4/12, 20 deg 1-2/12, 10 deg 0-1/12).
      Robust means and the whole QRS complex did not do better.
- [x] HRSP of the induced power (Cohen & Donner 2013) plotted by default ('hep_tf_power').
- [ ] A better sample dataset than sub-047 (user: "maybe we should find a better one"): screen
      other public EEG + ECG resting datasets (OpenNeuro) with the 2026-11-03 criteria (clear
      late HEP, flat pre-R, low frontal-pole activity, late map unlike the QRS map).
- [x] Surrogate reference (review 2026-11-04; user, 2026-11-05: keep the cycle average): the null compares
      each latency with the cardiac-cycle average, so a large locked deflection makes the other
      latencies differ too. Measured 2026-11-05: a pre-R reference fixes it in simulation
      (tests/sim_surrogate_reference.m: 2-9% vs 15-47% false positives) but no fixed pre-R
      window is quiet on real recordings (tests/check_surrogate_reference.m: P wave at -100 to
      -168 ms; sub-047 loses its late HEP with [-150 -50], sub-062 gets 44-308 ms with
      [-300 -200]). Recommendation: keep the cycle average (README explains).
- [x] Runtime of the family-wise HRSP/HRPC statistics (341 s vs 207 s, 2026-11-04): TFCE per
      height on the series that reach it, TF surrogates batched (2026-11-05, identical output):
      16 ch 170 -> 103 s; TFCE still 28% of it (74 s without), so the max statistic is the
      default for all measures (user, 2026-11-05) and TFCE is computed only when chosen.
- [ ] Line-noise detection: only 50-Hz recordings (66) and a synthetic 60-Hz one (the sample
      with 50 Hz notched out and 60 Hz added, run_tutorial_headless 'hep_ica_epochs_60hz');
      a real 60-Hz recording is still to find.
- [x] 'ecg_regression': the ICA copy keeps the cardiac field (ICLabel sees it, the data no
      longer have it). Harmless (tests/compare_regression_copy.m, 4 recordings, 2026-11-05):
      ROI HEP within 0.06 uV with the copy regressed too, heart-locked residual within .01 uV.
- [x] Tests for 'ica_epochs', 'ica_srate', 'ica_lowpass' > 0, 'ica_asr' 0 and the
      Picard-to-Infomax fallback (Picard stub with NaN weights): run_tutorial_headless
      'hep_ica_epochs_60hz', 'hep_ica_lowpass_noasr'. The copy/data mismatch fallback is a
      guard that no path reaches (same channels and epochs by construction): not tested.
- [x] Memory at 2048 Hz (user, 2026-11-05): no chunking; data above 1024 Hz must be
      downsampled (error on the command line, question in the GUI, 512 Hz recommended), above
      512 Hz a warning or question (run_checks).
- [x] 'rigid' surrogates removed (2026-11-05): still heartbeat-locked, they test condition
      contrasts, which BrainBeats does not do; 'hep_surrogate_mode' warns.
- [x] cfa_map offset rule (user, 2026-11-05): HEP without 'clean_eeg' stops on unfiltered EEG
      (run_checks: median |channel mean| > median 1-s SD; raw 105-1080, filtered < 0.6);
      cfa_map no longer removes offsets, it warns.
- [x] Paper on the CFA map vs posture and head orientation: dropped (user, 2026-11-05).

## Next release, after the current updates (user, 2026-11-05)
- Methods paper: the same public HEP dataset(s) (a well-known one or two, to be found)
  analyzed with different approaches, showing how each changes the results, and new
  guidelines: Alday's regressor at the group level, GLMs with inverse-variance weighting,
  covariate models (respiration, QRS amplitude, pre-R power), surrogate controls, IC vs scalp
  domain, building on the HEP_neurofeedback project; with the group-analysis tools released
  in BrainBeats if possible.
- README numbers with tolerances (tests/readme_numbers.m) once the program is updated.

## From the code review of v1.6 (fixed ones are in the README version history)
Scientific, may change results:
- Brain-heart coherence (experimental, out of the GUI): redesign at HRV time scales (NN series
  not high-passed at 0.5 Hz, MVAR at a matching rate); fit at the data rank with QR, not
  `(Y*Z')/(Z*Z')` and inv() in fdMVAR_5order, since cleaned EEG is rank-deficient.
  Redesign in progress on branch `brain-heart` (docs/brain-heart-plan.md there): MVAR at
  2 Hz on EEG band envelopes, NN and respiration; conditional GC; simulator first.
- [x] HRSP and the QRS (checked 2026-11-03): a heartbeat-locked broadband transient makes every
  point away from it a significant power decrease. Fixed 2026-11-04: the HRSP plotted is that
  of the induced power (HEP subtracted, Cohen & Donner 2013; 'hep_tf_power'): -0.08 dB and no
  FDR point instead of -1.9 dB and all (tests/test_compute_hep_tf.m check 8). Still open: the
  size of the residual on real data, where the QRS amplitude varies with breathing.
- [x] Surrogate beats can fall on epochs rejected by the cleaning: the rejected windows and the
  data outside the epoched range are passed to compute_hep_tf ('exclude'; user, 2026-11-05),
  surrogate trains refilled to the HEP's number of beats. Tutorial dataset null SD 0.43 vs 1.99
  uV (its first 2 s hold 5,000-uV artifacts); sub-047 +4%, sub-070 +1%.
- [x] HRV: LF and HF from the same 120-s windows, spread over the whole recording (2026-11-05,
  tests/compare_hrv_windows.m). 60/90-s windows let VLF into LF (1/f-like synthetic:
  Lomb-Scargle +24/+13%, welch +28/+7%; 120 s: +10/+2%); the 34-s HF windows let LF into HF
  (+14% vs +5%).
- Lomb-Scargle has no taper: VLF power adds to LF (+7-10% on synthetic series with VLF, welch
  +/-2%). Try detrending or tapering each window; compare with tests/compare_hrv_windows.m.
- [x] clean_rr ectopic rule: 'adaptive' (Lipponen & Tarvainen, 2019) is the default since
  2026-11-05 (user asked for it; tests/compare_ectopic_rules.m, 40 recordings): 92-96% vs 13-15%
  of injected premature beats at 0.80 x, RMSSD error 1-7 vs 13-23 ms; replaced 1 sinus and 2
  unclear beats in 12,800 (20% rule: none). 'ectopic_rule','ratio' keeps the v1.6 rule.
- [x] Ectopic rule on the MIT-BIH Arrhythmia Database (2026-11-06, tests/validate_ectopic_mitbih.m):
  0.66 false replacements per 1000 sinus beats; 76% of sparse premature beats replaced, 5-10%
  where over 15% of beats are abnormal. New warning when the cleaned series still alternates
  (RR_ALTERNATION > 10%): 14/44 records (all with AF/flutter or frequent ectopy), 0/65 healthy.
  Not caught: 207 (ventricular flutter, 153 false detections, RMSSD -162 ms), 232 (frequent
  atrial premature beats with long pauses).
- PAT: show the method in the HEP figure titles; check the censoring of the pairing in 'peaks'
  mode. (The 250 ms fallback is gone: without an ECG, 'ppg_transit' must be given or 'eeg'.)
- Validation: every number in the README comes from one subject. Run other ds003838 subjects,
  an annotated ECG/PPG benchmark and an HRV cross-check; v1.4 vs v1.6 table of changed outputs.
Robustness:
- [x] Absolute thresholds (2026-11-05): rank relative (data_rank: (n eps_single)^2 x the
  largest eigenvalue; same ranks on 4 recordings in uV, right in volts where 1e-7 gave 0-1);
  ECG flatline = no variation around the offset (the old check, after the band-pass, let a
  constant with an offset through); signal_filt index rounded for non-integer fs; PPG
  clipping relative to each beat's range per second (same 'Q' share on our data, 0.3-0.8%,
  in any unit; the old one went to 4-23% at x1e-3).
- [x] Fewer than 30 channels: no re-reference, with a warning that says so and suggests choosing
  one (user, 2026-11-05: leave it to the users).
- [x] Bad-epoch HF metric (2026-11-05, tests/compare_badepoch_metrics.m): analysis window only;
  HF = 90th percentile over channels of the RMS of the >45 Hz residual of the broadband ICA
  copy: 15/15 added muscle bursts on 5 of 6 recordings (old: 0-2/15).
- IBIs are computed after synthetic beats are dropped: compute them before.
- 'rr' mode cannot take a PAT: accept a numeric 'ppg_transit'; tests/test_rr_vs_ppg.m locks HEPs
  to pulses.
- Nonlinear EEG features are O(N^2) on whole recordings and on by default: use 30-60 s windows.
- Generic helper names on the path: package them in the redesign.
- Tests: most headless cases only check that the run finishes: add tests/readme_numbers.m (the
  numbers quoted in the README, with tolerances) and a synthetic HRV and EEG features test.
Documentation:
- References checked 2026-11-03: Park et al. (2018, Cerebral Cortex 28:2351, surrogate R-peaks
  shifted by -500 to +500 ms, after Park et al., 2014) and Lee et al. (2024, NeuroImage 299:120797)
  are right; Ros et al. (2025) is a bioRxiv preprint (GEDAI); Steinfath et al. (2026), the
  systematic review of heartbeat-evoked responses (bioRxiv 2025.08.08.668923), check the journal
  version before the release; the pseudotrial paper (Imaging Neuroscience, doi:10.1162/IMAG.a.30):
  add its authors to the README.

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

Decided (user, 2026-11-06), for v2.0 after the v1.7 release, on its own branch:
- GUI: successive windows (main window: data and heart signal; then preprocessing; then
  measures with Options buttons), with Back.
- Outputs: [EEG, com, HEP]: the continuous cleaned dataset and the HEP epochs; in EEGLAB the
  HEP becomes a new dataset. Legacy 'analysis','hep' keeps returning the epochs first.
- Heart-component removal reuses the cleaning ICA; off by default for every measure ('ica':
  ICLabel heart >= .75; the features path's fixed .99 rule goes).
- One shared cleaning for all measures (one ICA on the continuous 1-Hz copy), then a second
  pass by data type: bad epochs for the HEP; ASR in reconstruct mode after the ICA for the
  continuous data (cutoff to validate: 10/15/20/30 SD, user suggested 15).
- A measure's results must not depend on which other measures are selected.

## Next (other)
- Heart artifact removal ('heart_removal'): default 'none' (user, 2026-11-02: "add all of them and
  make it clear ... we actually recommend NOT to do them"); options 'ica', 'ica_caracas', 'ica_ecg',
  'ecg_regression', 'ssp'. 6 ds003838 subjects (tests/compare_heart_removal.m, synthetic on): the
  three IC-selection methods found no heart IC on any subject (ICLabel; CARACAS; |r| with the ECG
  below .6, .55 on the sample); regression removed 29% (3-70) of the QRS amplitude and kept 80% of a
  synthetic response; SSP 52% (27-73) and 57% (39-76). Open:
  - CARACAS on the sample: the most ECG-like IC fails only skewness (0.90 vs 1). Longer or
    higher-CFA recordings (e.g. more channels, supine) would test the IC methods better.
  - A realistic synthetic HEP (variable amplitude and latency, cortical topography).
  - Report the ECG's heartbeat-locked average per condition next to the HEP (README advises it).
- Surrogates also for HEP contrasts between conditions (e.g. per-condition beat subsets).
- 'individualized' EEG bands: only alpha is individualized (per-channel peaks are unreliable
  for the other bands).
- get_RR fails on recordings shorter than ~3 s (filter length).
- PAT without an ECG (tests/compare_ppg_transit.m, 6 ds003838 subjects, 2026-11-02): the EEG
  estimate was wrong on 2/6 (628 and 116 ms vs 408 and 392 ms from the ECG) with ratios 2.2/2.1,
  like a correct one (2.1); the 250 ms default was wrong on 6/6 (true 384-424 ms). With an ECG the
  correction works (r .52-.86 with the ECG-locked HEP). Decided (2026-11-02): no silent correction
  without an ECG ('auto' errors with the EEG estimate as a suggestion; 'eeg' only on request).
  After the PPG polarity fix (2026-11-03, pulse onsets): PAT 92-112 ms, ECG correction r .68-.92;
  EEG estimate right on 3, wrong on 1 (404 vs 104 ms, ratio 2.05), not clear enough on 2.
- PPG pulse onset: detected on the 0.5-3 Hz signal (R + 92-112 ms after the polarity fix); locate
  it on a less smoothed signal (e.g. maximum of the second derivative) and compare the HEP match.
- Python port (../BrainBeats-py): carry over the 2026-11-03 changes (PPG polarity, ectopic rule,
  Lomb-Scargle ms^2 default and grid, half-open EEG bands, asymmetry sign, qeeg paf, no 'mode')
  and the v1.7 changes (HRV LF/HF from the same 120-s windows, bad-epoch metric, relative
  thresholds, surrogate exclusions, max statistic default, RR_ALTERNATION warning).
- EEG aperiodic activity (user, 2026-11-03): add the aperiodic exponent and offset, and the
  spectrum and band power after removing the aperiodic component compared with before, as in
  ASCENT (github.com/amisepa/Ascent, specparam; offer to install the plugin on the fly when
  missing, as for GEDAI).
- 'keep_heart' stores the raw heart channel (sample ECG has a large DC offset): filter it?
