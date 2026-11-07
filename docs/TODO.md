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
- [ ] ICA quality for HEP (user: improve the decomposition rather than lower the ICLabel eye
      threshold, .90 kept): compare Picard (used for the README figures and the screening) with
      extended Infomax (the default, 'icamethod' 2), with and without ASR on the ICA copy (item
      below). Measures: ICLabel brain/eye IC counts, eye activity left in the kept epochs,
      ROI HEP and its surrogate test.
- [ ] Surrogate statistics more conservative than per-latency FDR (user): max-statistic (tmax,
      FWER: null = max |z| over latencies of each surrogate) and 1-D TFCE over time, null from
      the surrogates' TFCE maxima. Needs more surrogates (100 gives a coarse 95th percentile of
      maxima; ~1000 for the HEP, cheap without HRSP/HRPC). Validate on noise (FWER ~5%) and on
      sub-047 before switching the default.
- [ ] IC figure: fewer plots. Clustering within one recording would merge independent sources;
      rather sort by heartbeat locking (max |z|) and optionally show only ICs that pass the
      surrogate test or are classified brain. Across subjects, cluster ICs in an EEGLAB STUDY
      (scalp maps + HEP).
- [ ] README: hep_ecg and hep_ppg look the same by design (the PPG beats are shifted back by the
      pulse arrival time measured from the ECG). A single overlay of the ROI HEP from the ECG,
      the corrected PPG and the uncorrected PPG would show that instead.
- [ ] plot_NN labels the ECG in uV; ds003838 stores the ECG in raw amplifier units (50,000-
      730,000 at the 99.5th percentile across 21 recordings, the tutorial dataset included).
      Label from EEG.chanlocs/units when known, else 'a.u.'.
- [ ] More ds003838 subjects if a clearer HEP is needed: 22 of 65 usable screened (sample, 033-036,
      038-054); 055-098 not yet.
- [ ] HEP cleaning, ASR on the ICA copy (user, 2026-11-03): ICA is already fitted on a 1-Hz
      high-passed copy (filtered continuous, then epoched; bad epochs removed from both) and the
      weights applied to the 0.5-Hz data, but the HEP path has no ASR (bad epochs only). Try ASR
      on the 1-Hz continuous copy before epoching it, for the ICA fit only (the HEP data are not
      modified by ASR), with a very lax cutoff (~100 SD, user) so ocular activity stays in the
      copy and ICA can model it; compare a few cutoffs (e.g. 50, 100). Measure first on the 22 screened ds003838 recordings: ICLabel brain/eye
      IC counts, frontal-pole activity left in the kept epochs, ROI HEP and its surrogate test,
      with and without. Same rank question: ICA's rank is counted with an absolute 1e-7
      threshold on the covariance eigenvalues (uV^2), in clean_eeg and run_HEP.

## From the code review of v1.6 (fixed ones are in the README version history)
Scientific, may change results:
- [x] Brain-heart coherence redesigned as 'coupling' (branch brain-heart, docs/brain-heart-plan.md):
  beat-to-beat series, VARX with respiration exogenous, state-space GC, tested on simulations.
  Done: ds003838 cohort (61 usable of 64, all variants), ECG-derived respiration vs a belt
  (38 participants), cardiac field options on real data ('evoked' kept), calibration at 200
  recordings per scenario with fresh seeds, segment-swap null on real data, report figure.
  Numbers in README section 3. Release plan (user, 2026-11-06): v1.7 = main without coupling;
  v1.8 = this branch merged after the v1.7 release; v2.0 = the step-based redesign (its own
  branch). For the v1.8 merge:
  - [ ] Merge main (v1.7) into this branch here and resolve: brainbeats_process.m,
        getparams_cmd.m, clean_eeg.m, tests/make_readme_figures.m, docs/TODO.md (trial
        merge, 2026-11-06).
  - [ ] v1.7 changes the heart series (adaptive ectopic rule by default) and adds the
        rr_alternation warning (frequent ectopic beats, atrial fibrillation): carry that
        warning into the coupling report (no sinus rhythm, no heart-period GC), then rerun
        test_coupling_engine/sim/stats and the ds003838 cohort on the merged code.
  - [ ] Label transfer entropy, time-reversed GC, DC/GPDC as exploratory in the README (TE
        rejected 10% heart -> alpha in simulation; the report uses the F-test and q only).
  Still open:
  - [x] F-tests ~1-2 points liberal on real nulls: now HC3-robust with the degrees of freedom
        less the detrended runs (5.8% on both real nulls, was 6.1-7.2%; README section 3).
  - [ ] Still ~1 point above 5% in places: an FDR-significant link in 11% of the ds003838
        null models (4% in the belt dataset; 5% expected under the complete null); 6.5% on
        HEP-NF swaps with each segment cut in two. Watch it in the paper's null; no cause
        found beyond heteroskedasticity.
  - [ ] Respiratory phase inputs ('resp_phase' in get_coupling_series): helped EEG -> heart
        calibration in simulation only; test on the belt dataset before exposing it.
  - [x] Heart-period equation failed Ljung-Box (p < .05) in 58% (belt dataset) and 40%
        (ds003838) of the recordings: own lags are now added until the residuals are white
        (0-1% fail; nulls 5.4% and 5.7%, real pairings unchanged). Log NN, HR, longer BIC
        caps, cubic detrending and longer respiration lags did not help.
  - [ ] Belt-dataset heart residuals are heavy-tailed (kurtosis ~14 vs ~4 in ds003838):
        residual artifacts in its RR series (cached before the v1.7 adaptive ectopic rule)?
        Check after the v1.8 merge.
  - [ ] Interpolate single missed PPG pulses: scattered misses break the PPG series into
        pieces < 20 beats, and the PPG proxies are then left out (3 of 61 in ds003838).
  - [ ] A dataset with a belt and long recordings for the paper (MIT-BIH Polysomnographic,
        CAP sleep); a paced-breathing recording.
  - [ ] Two runs of the README coupling figure (sub-047, same settings) differed at the
        margin (frontal theta -> heart mediated in one, uncorrected only in the other; CFA
        max 4.5% vs 4.4%): the cleaning is not deterministic. Find which step, and seed it.
  - [ ] GUI entry, Python port (MNE-Connectivity gc/gc_tr, Tigramite PCMCI+).
- HRSP and the QRS (checked 2026-11-03): a heartbeat-locked broadband transient 5 x the noise SD
  makes every point at 350-600 ms a significant power decrease (-1.7 dB) against the 'random'
  surrogates (none at 2 x); normalizing by the continuous data gives the same (-1.9 dB), since it
  contains every QRS. Possible fix: a QRS-free null (surrogate power from samples away from the
  R-peaks); impossible at low frequencies (a 5-cycle wavelet at 4 Hz spans a whole cycle). Size
  of the effect on real data: measure on ds003838 (HRSP at 350-600 ms vs the CFA amplitude).
- Surrogate beats can fall on epochs rejected by the cleaning (wider null, depends on the
  artifact load): pass the rejected spans to compute_hep_tf as exclusions.
- HRV: LF (125-s windows) and HF (34-s windows) come from windows of different lengths; with ms^2
  (default since 2026-11-03) both estimate the same spectrum, but LF and HF from one window of
  2 min or more would be simpler. LF/(LF+HF) = 'hrv_norm' when VLF is not estimated.
- clean_rr ectopic rule (20%): detected 1-6 of 6 injected premature beats at 0.70 x per
  recording, 1-4 at 0.75 x; looser rules flagged beats in the high-RMSSD sub-033. An adaptive
  rule (Lipponen & Tarvainen, 2019) could do better: compare on an annotated arrhythmia database.
- PAT: show the method in the HEP figure titles; check the censoring of the pairing in 'peaks'
  mode. (The 250 ms fallback is gone: without an ECG, 'ppg_transit' must be given or 'eeg'.)
- Validation: every number in the README comes from one subject. Run other ds003838 subjects,
  an annotated ECG/PPG benchmark and an HRV cross-check; v1.4 vs v1.6 table of changed outputs.
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

Decisions needed: tabs vs successive windows; second output for the HEP dataset; whether
heart-component removal keeps the heart channel in the ICA (as 'rm_heart' did) or
reuses the cleaning ICA (faster, one decomposition).

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
  Lomb-Scargle ms^2 default and grid, half-open EEG bands, asymmetry sign, qeeg paf, no 'mode').
- EEG aperiodic activity (user, 2026-11-03): add the aperiodic exponent and offset, and the
  spectrum and band power after removing the aperiodic component compared with before, as in
  ASCENT (github.com/amisepa/Ascent, specparam; offer to install the plugin on the fly when
  missing, as for GEDAI).
- 'keep_heart' stores the raw heart channel (sample ECG has a large DC offset): filter it?
