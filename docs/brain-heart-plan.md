# Brain-heart coupling: redesign plan (branch `brain-heart`)

Replaces the experimental 'coherence' mode. Merge into main only after the simulation
tiers below pass. Every default is provisional until validated.

## Why
The current mode fits an order-5 MVAR to raw EEG at the EEG rate plus the NN series
band-passed 0.5-40 Hz: HRV (0.04-0.4 Hz) is filtered out and 5 lags span 5-20 ms, so it
models the ECG trace, not autonomic coupling.

## Questions the module answers
| Question | Primary | Check |
|---|---|---|
| Coupling, how strong | COH; spectral GC (bits) over LF and HF | Gaussian-copula MI |
| Direction | conditional GC, both directions, own test each | time-reversed GC; GPDC/DC spectra |
| Direct or via respiration / other EEG | conditional vs pairwise GC; pCOH vs COH | CFTd vs CFTf surrogates |
| Instantaneous | zero-lag residual partial correlation | |
| Nonlinear beyond linear | GCMI conditional TE vs linear TE (= GC/2) | KSG TE excess (optional, JIDT) |

DC and DTF are total measures; PDC, GPDC, pCOH and conditional GC are direct. PDC/DC are
normalized ratios: not comparable across models with different variable sets. GC in bits is.

## Phases
0. **Simulator first** (tests/): tier 1 linear VAR with analytic GC/DC/PDC (null, one-way,
   two-way, chain via respiration, respiration driving both, hidden driver); tier 3 end to
   end: EEG envelope + 1/f, IPFM heartbeat generator with set gains and lags, respiration
   driving vagal input and EEG, synthetic ECG, QRS-locked CFA injected into EEG, ectopic
   beats. Scoring: Type-I error (also under a common driver), power/ROC vs length (2, 5,
   10 min), direction confusion matrix, effect-size bias, null calibration.
1. **Series**: EEG band envelopes (band-pass, Hilbert, anti-aliased decimation) for 2-4
   regions or components, after CFA removal, average reference; NN interpolated, detrended,
   no high-pass; respiration from a channel, else derived from ECG/PPG and flagged. Common
   2-Hz grid (4 Hz optional). 5-min stationary windows by default.
2. **Linear engine**: VAR, BIC order over lags up to 10 s, rank-safe fit (QR), samples per
   coefficient >= 10 or LASSO/block GC; state-space conditional and pairwise GC (time and
   spectral); COH, pCOH, GPDC, DC from the same fit; residual whiteness and zero-lag residual
   correlation reported. MVGC2 (GPL-3.0) as an optional engine if mature, else own code.
3. **Statistics**: F/chi2 for time-domain GC; causal Fourier surrogates (CFTd direct, CFTf
   total) for spectral measures; circular-shift or IAAFT for MI-type measures; never
   independent shuffles. Group level: max-statistic or FDR across region x band x direction.
4. **Nonlinear complement**: GCMI conditional TE on the same conditioning set (gcmi,
   GPL-3.0); optional KSG conditional TE via JIDT with linear-preserving surrogates.
5. **Outputs**: one row per source, target, band, measure: value, p, null type,
   conditioning set, order, samples/coef, whiteness, zero-lag r, TRGC agreement,
   respiration source (measured/derived/absent), change without CFA removal. Two-layer
   graph (direct vs mediated); spectra on the 0-1 Hz modulation axis. Warnings when
   respiration is absent, windows are short, or TRGC disagrees.
6. **Real data**: MIT-BIH Polysomnographic (EEG, ECG, respiration, BP), CAP sleep, ds003838
   (no respiration), VitalDB (anaesthesia); negative controls for direction: arousals,
   ictal tachycardia (Siena). A paced vs spontaneous breathing recording would close the
   gap no open dataset covers.

## Decisions from measurements (2026-11-05)
- Respiration is an exogenous input (VARX, lags 0..p), not a VAR variable: on the tier-3
  simulator a smooth breathing signal needed univariate order 28-30 with residual
  Ljung-Box p ~1e-175 and pushed the whole model to the order cap (9 at 2 Hz) with
  spectral radius 0.997. As an input: order 4, radius ~0.9, white EEG residuals, and the
  planted resp -> heart / resp -> alpha links recovered. GC among EEG and heart given a
  known input follows from A and V alone.
- One BIC order per equation (Hsiao, 1981): EEG envelopes need 1-2 lags, the interpolated
  NN series ~6 at 2 Hz. Heart residuals still not white (Ljung-Box p 1e-6 to 1e-12):
  watch the calibration; own-lag vs cross-lag orders are the next refinement if needed.
- Own and cross lag orders separately per equation (var_design / var_order crit.P): on the
  2-Hz grid the NN equation took 8-11 own lags and 1 cross lag, but its residuals still
  failed Ljung-Box in 88% of null simulations and alpha->heart / heart->beta F-tests
  rejected 10-15%: a spline-interpolated NN series sampled faster than the heart rate
  carries between-beat structure an AR model cannot whiten.
- Beat-to-beat grid by default (as in cardiovascular variability analysis): one sample
  per heartbeat, NN itself, mean band power and respiration over each RR interval. Ljung-Box
  failures 5-7.5% in every equation (nominal), NN own order ~5. Frequencies in Hz at the
  mean heart period. A uniform grid ('coupling_grid','uniform') is kept as an option.
- Pairwise = a separate bivariate VAR of the band and the heart (no respiration, no other
  bands); total coherence from that model, partial coherence from the full model.
- Cardiac field: 'evoked' stays the default. In simulation (13 channels) 'both' (SSP +
  evoked) calibrated best (3-6% vs heart -> beta 8-9%) but cut power (b2h 92%, h2b 70% vs
  100%, 98%); in ds003838 (63 channels) the options gave the same rates (heart -> beta 8.6%
  evoked, 8.2% none, 14.3% both).
- Detrending per run of complete samples: one detrend over pooled 90-s segments let block
  offsets look like lagged prediction. The 25% significant tests that remained in the belt
  dataset survived a segment-swap null (7.7%): real within-segment co-fluctuation.
- PPG: respiration from the PPG is as good as from the ECG (HF coherence with NN 0.86 /
  0.73 vs 0.87) and all three together raise resp -> heart from 94% to 98% without
  changing EEG <-> heart; pulse arrival time and amplitude capture the cardiovascular loop
  but no EEG link above chance at rest: options, off by default. A PPG series that would
  cut the usable data by more than half is dropped with a warning.
- FDR (q) on the F-test p: the bootstrap's resolution, 1/(nboot+1), is too coarse for ~32
  tests. A model that cannot be fitted in any region is an error (it used to crash later).
- F-tests robust to heteroskedasticity (HC3 Wald) with the residual degrees of freedom less
  2 per detrended run. Classical F on real nulls (EEG of one segment or recording with the
  heart of another): 6.1% (HEP-NF), 7.2% (ds003838); robust: 5.8% and 5.8%; simulation
  5.4 -> 5.3%; power unchanged. Not the cause (measured): order selection (cross order 1
  in 99% of tests), residual autocorrelation (Bartlett/Dutilleul inflation 1.00), BIC
  underfitting, onset-aligned trends. A wild bootstrap did no better than HC3.

Comparators for the paper only, not shipped: MIC at 1 Hz and lag 0, SDG (Catrambone code),
PCMI (order <= 3, circular-shift null per direction), pairwise GC.

Python (Tigramite PCMCI+, IDTxl, MNE-Connectivity gc/gc_tr) goes to the python-port branch.
