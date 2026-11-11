# BrainBeats TODO

Open items only. What was done, with its numbers, is in the README version history and in the
headers of the test and validation scripts.

## Before the v1.7 release
- [ ] Update the installed plugin first: eeglab/plugins/BrainBeats (a clone, what EEGLAB
      runs) was at 9c34514 with 26 local changes on 2026-11-11, so the GUI would run older code.
- [ ] Click through the GUI in the MATLAB desktop: 'BrainBeats' in the EEGLAB menu bar, main
      window (Load dataset / Sample dataset, analysis, heart channels and the '...' list,
      'Perform analysis on'), parameters window for HEP and for features (heart artifact
      removal options, PPG pulse arrival time options), and that the HEP output appears as a
      new dataset. Automated tests check the parameters and the history command but draw no
      window.
- [ ] Look at the README figures in the MATLAB desktop once (made by
      `tests/make_readme_figures`).
- [ ] ASCENT 2.1 or later on the EEGLAB extension manager, under the name 'ASCENT': BrainBeats
      installs it from there for the fuzzy entropy and fractal dimension (it is only on GitHub
      now). Rename ASCENT's compute_psd first: BrainBeats 1.5-1.6 (1.5 is on the extension
      manager) has a compute_psd with another signature, and whichever comes first on the path
      breaks the other plugin.
- [ ] Release: date in the README version history, tag `eeglab_brainbeats_v1.7`, GitHub
      release (zip), EEGLAB extension manager.

## v1.8: brain-heart coupling
- Merge branch `brain-heart` (coupling replaces the experimental coherence analysis; plan in
  docs/brain-heart-plan.md on that branch) after the v1.7 review fixes are committed, then the
  release check.
- Conflicts to expect: brainbeats_process.m (help: main dropped 'asr_mem' and documents the
  ASCENT features; the coupling section and history), clean_eeg.m (the branch renames
  'coherence' to 'coupling' in the ASR branch; main replaced the RAM-sized ASR call there),
  getparams_cmd.m ('asr_mem' now warns), README.md, docs/TODO.md, brainbeats_tutorial.m,
  tests/run_tutorial_headless.m, tests/make_readme_figures.m. The branch calls none of the
  renamed or removed functions (compute_psd, scrollplot, compute_fe, fractal_volatility).
- The coupling cleaning gets main's fixed ASR chunks: the branch's cohort numbers and
  ds003838/coupling_cache were made with RAM-sized chunks, which change the cleaned data
  between runs. Rebuild the cache and rerun the cohort before quoting numbers.

## v2.0: step-based redesign (decided 2026-11-06, own branch, after v1.8)
- Steps instead of exclusive analyses: heart (beats, QA, RR cleaning, PPG transit), EEG
  cleaning, optional heart-component removal, then any subset of measures (HEP with
  HRSP/HRPC and surrogates, HRV features, EEG features, coupling) in one run. 'measures' key;
  the 'analysis' values map onto it so existing scripts keep working.
- A measure's results must not depend on which other measures are selected.
- One shared cleaning: one ICA, fitted on the continuous 1-Hz copy (lax ASR, reconstruct
  mode), the same components removed for every measure; then a second pass by data type: bad
  epochs for the HEP, ASR in reconstruct mode after the ICA for the continuous data (cutoff to
  validate on 10/15/20/30 SD; 15 suggested). Low-pass per output (continuous 40 Hz, HEP 30 Hz).
- Heart-component removal reuses the cleaning ICA; off by default for every measure ('ica':
  ICLabel heart >= .75; the features path's fixed .99 heart rule goes).
- Outputs: [EEG, com, HEP]: the continuous cleaned dataset and the HEP epochs (a new dataset in
  EEGLAB). 'analysis','hep' keeps returning the epochs first.
- GUI: successive windows (data and heart signal; preprocessing; measures with an Options
  button each), with Back.
- Colormap: the red-blue diverging map ('RdBu_r') for every signed map (HEP topographies and
  images, QRS map, z-scores, HRSP in dB), limits symmetric around 0; a sequential map for
  unsigned ones (HRPC, power, probabilities). One helper for every figure.
- Validation before switching: HEP-only v1.7 vs v2.0 (expected identical), EEG features v1.7 vs
  v2.0 (report the change), on the two sample datasets and ds003838 sub-033/040/062/070.
- Package the generic helper names (filter_fast, plot_topo, ...) to avoid path clashes.

## Next release, after v2.0
- Methods paper: one or two well-known public HEP datasets analyzed with different approaches,
  showing how each changes the results, with guidelines: Alday's regressor at the group level,
  GLMs with inverse-variance weighting, covariates (respiration, QRS amplitude, pre-R power),
  surrogate controls, IC vs scalp domain; group-analysis tools released in BrainBeats.
- `tests/readme_numbers.m`: the numbers quoted in the README, with tolerances.

## Science and validation
- IC figure: fewer plots, sorted by heartbeat locking (max |z|), optionally only the ICs that
  pass the surrogate test or are classified brain; across subjects, cluster ICs in a STUDY.
- A better sample dataset than sub-047: screen other public EEG + ECG resting datasets (clear
  late HEP, flat pre-R, low frontal-pole activity, late map unlike the QRS map).
- Line-noise detection: tested on 50-Hz recordings and a synthetic 60-Hz one only; find a real
  60-Hz recording.
- HRSP: size of the induced-power residual on real data, where the QRS amplitude varies with
  breathing.
- Lomb-Scargle has no taper: VLF power adds to LF (+7-10% on synthetic series with VLF, welch
  +/-2%). Try detrending or tapering each window (tests/compare_hrv_windows.m).
- Ectopic rule: frequent ectopy is only flagged (RR_ALTERNATION); MIT-BIH records 207
  (ventricular flutter, false detections) and 232 (atrial premature beats with long pauses)
  are not.
- Heart artifact removal: test the IC methods on recordings with a larger cardiac field (more
  channels, supine); a realistic synthetic HEP (variable amplitude and latency, cortical
  topography); report the ECG's heartbeat-locked average per condition next to the HEP.
- Surrogates for HEP contrasts between conditions (per-condition beat subsets).
- PAT: show the method in the HEP figure titles; check the pairing in 'peaks' mode; locate the
  pulse onset on a less smoothed signal (e.g. maximum of the second derivative).
- Validation beyond the sample data: other ds003838 subjects, an annotated ECG/PPG benchmark,
  an HRV cross-check, and a table of outputs changed since v1.4.
- 'individualized' EEG bands: only alpha is individualized.
- EEG aperiodic activity: exponent and offset, and the spectrum and band power with the
  aperiodic component removed, as in ASCENT (specparam; offer to install it when missing).

## Robustness
- IBIs are computed after synthetic beats are dropped: compute them before.
- 'rr' mode cannot take a PAT: accept a numeric 'ppg_transit'.
- Nonlinear EEG features are O(N^2) on whole recordings and on by default (ASCENT's fuzzy
  entropy is twice as fast as v1.6's, still O(N^2)): use 30-60 s windows.

## Python port (branch python-port)
- Carry over the v1.7 changes: PPG polarity, adaptive ectopic rule and RR_ALTERNATION,
  Lomb-Scargle in ms^2 and its grid, half-open EEG bands, asymmetry sign, qEEG paf, no 'mode',
  HRV LF/HF from the same 120-s windows, bad-epoch metric, relative thresholds, surrogate
  exclusions, max statistic default, pNN50 > 50 ms on differences rounded to 1 us, median
  and SEF90 from the bin edges, 'hrv_norm' NaN without LF, sub-sample beat times (vertex
  of the parabola through each peak or onset and its neighbors) for RR and RR_t, PPG
  polarity from the quantile skewness of the slope (5th/95th percentiles, warning below 0.2),
  'keep_heart' channel filtered like the EEG, get_RR minimum durations (ECG 3 s, PPG 12 s),
  cleaning plots in 'a.u.'.

## Decisions (with the reason, so they are not reopened)
- 2026-11-02: heart artifact not removed by default; ICA, CARACAS, ECG correlation, ECG
  regression and SSP are options to compare (each removed part of a synthetic response).
- 2026-11-02: no silent PPG pulse-arrival-time correction without an ECG (the EEG estimate was
  wrong on 1-2 of 6 subjects).
- 2026-11-03: ICA for HEP on a continuous 1-Hz copy with lax ASR (80); ICLabel eye threshold .90
  kept.
- 2026-11-04: no recording-level baseline regression (a centered regressor cannot change the
  average HEP); pre-R values exported for group models instead.
- 2026-11-05: surrogate reference is the cardiac-cycle average (no pre-R window is quiet on real
  recordings); max statistic by default (TFCE costs 28% more time).
- 2026-11-05: fewer than 30 channels: no re-reference, a warning.
- 2026-11-05: posture/head-orientation paper on the CFA map dropped.
