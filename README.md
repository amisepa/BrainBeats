<!-- <p align="center"> -->
# BrainBeats
<!-- </p> -->

<p align="center" width="100%">
  <img width="30%" alt="BrainBeats logo"
       src="https://raw.githubusercontent.com/amisepa/BrainBeats/main/brainbeats_logo2.png">
</p>


The BrainBeats toolbox, implemented as an EEGLAB plugin, allows joint processing and analysis of EEG and cardiovascular signals (ECG and PPG) for brain-heart interplay research. Both the general user interface (GUI) and command line are supported (see tutorial). BrainBeats currently supports: 1) Heartbeat-evoked potentials (HEP) and oscillations (HEO); 2) Extraction of EEG and HRV features; 3) Extraction of heart artifacts from EEG signals; 4) brain-heart coherence.

 
## 4 METHODS AVAILABLE

<p align="center" width="100%">
  <img width="50%" alt="BrainBeats diagram"
       src="https://raw.githubusercontent.com/amisepa/BrainBeats/main/figures/diagram.png">
</p>


1) Process EEG data for heartbeat-evoked potentials (HEP) analysis using ECG or PPG signals. Steps include signal processing of EEG and cardiovascular signals, inserting R-peak markers into the EEG data, segmentation around the R-peaks (-300 to 600 ms by default, or a window adapted to the subject's heart rate), rejection of heartbeats whose next QRS would fall in the epoch and of outlier inter-beat intervals, optional regression-based baseline correction (Alday, 2019), and heartbeat-related spectral perturbations (HRSP). With PPG, the pulse arrival time can be corrected, from an ECG channel or a known delay. Beat times detected with another tool can be provided directly ('heart_signal','rr' with 'beat_latencies').


<p align="center">
    Example of HEP at the subject level, obtained from simultaneous EEG-ECG signals (the cardiac field artifact was preserved here for illustration).
</p>
<p align="center" width="100%">
  <img width="50%" alt="BrainBeats fig11"
       src="https://raw.githubusercontent.com/amisepa/BrainBeats/main/figures/fig11.png">
</p>

<p align="center">
    Example of HEP at the subject level, obtained from simultaneous EEG-PPG signals (note that with PPG, we must correct for the delay between the electrical and mechanical cardiac events so that the estimated heartbeat times correspond to the R-peaks of an ECG; ~200-400 ms depending on PPG sensor location, the subject, and if you use the onset or peak of the pulse wave as the marker).
</p>
<p align="center" width="100%">
  <img width="50%" alt="BrainBeats fig17"
       src="https://raw.githubusercontent.com/amisepa/BrainBeats/main/figures/fig17.png">
</p>

2) Extract EEG and HRV features from continuous data in the time, frequency, and nonlinear domains. 
    - HRV time domain: SDNN, RMSSD, pNN50.
    - HRV frequency domain: VLF-power, ULF-power, LF-power, HF-power, LF:HF ratio, Total power. 
    - HRV nonlinear domain: Poincare, fuzzy entropy, fractal dimension, PRSA. 
    
    - EEG frequency domain: average band power (delta, theta, alpha, beta, gamma), individual alpha frequency (IAF), alpha asymmetry.
    - EEG nonlinear domain: fuzzy entropy, fractal dimension


<p align="center">
    Example of power spectral density (PSD) estimated from HRV and EEG data
</p>
<p align="center" width="100%">
  <img width="50%" alt="BrainBeats fig21"
       src="https://raw.githubusercontent.com/amisepa/BrainBeats/main/figures/fig21.png">
</p>

<p align="center">
    Example of EEG features extracted from sample dataset
</p>
<p align="center" width="100%">
  <img width="50%" alt="BrainBeats fig22"
       src="https://raw.githubusercontent.com/amisepa/BrainBeats/main/figures/fig22.png">
</p>

3) Remove heart components from EEG signals using ICA and ICLabel.
   
<p align="center">
    Example of extraction of cardiovascular components from EEG signals
</p>
<p align="center" width="100%">
  <img width="50%" alt="BrainBeats fig27"
       src="https://raw.githubusercontent.com/amisepa/BrainBeats/main/figures/fig27.png">
</p>

4) Compute brain-heart coherence (beta version, please test and give feedback)
   
<p align="center">
    Example of several brain-heart coherence measures computed with BrainBeats from simultaneous EEG and ECG signals
</p>
<p align="center" width="100%">
  <img width="50%" alt="BrainBeats coh_all"
       src="https://raw.githubusercontent.com/amisepa/BrainBeats/main/figures/coherence_allfreqs.png">
</p>

<p align="center">
    Scalp topography showing scalp regions coherent with ECG signal for each frequency band
</p>
<p align="center" width="100%">
  <img width="50%" alt="BrainBeats coh_topo"
       src="https://raw.githubusercontent.com/amisepa/BrainBeats/main/figures/coherence_topo.png">
</p>


## Requirements

- MATLAB installed (https://www.mathworks.com/downloads), with the Signal Processing and the Statistics and Machine Learning toolboxes. The Parallel Computing toolbox is optional ('parpool' option).
- EEGLAB installed (https://github.com/sccn/eeglab), with the clean_rawdata, ICLabel and firfilt plugins (included in EEGLAB by default). PICARD (fast ICA) and REST (infinity reference) are installed automatically when these options are selected.
- Some data containing EEG and cardiovascular signals (ECG or PPG) within the same file (i.e. recorded simultaneously).
  Or use the tutorial dataset provided in this repository located in the "sample_data" folder. Source: sub-32 in https://nemar.org/dataexplorer/detail?dataset_id=ds003838

## Step-by-step tutorial

JoVe video tutorial: https://www.jove.com/v/65829/author-spotlight-advancing-study-brain-heart-interplay-with

JoVe publication tutorial: https://www.jove.com/t/65829/brainbeats-as-an-open-source-eeglab-plugin-to-jointly-analyze-eeg


## Version history

v1.6 (9/2026) - Many fixes and a few additions (see below). Results of some features change:
- Heartbeats: improved R-peak detection and correction (get_RR, clean_rr: a removed false beat now merges its interval with the next one; missing beats are inserted in gaps), per-electrode outputs with several heart channels, new R-peak QA report (rpeak_qa: flags beats detected on the wrong wave, e.g. T-wave), ECG artifact detector (detect_ecg_artifacts), PPG 'ppg_detect_mode' ('valleys' or 'peaks').
- New regression-based baseline correction of HEP epochs ('hep_baseline','regression'; Alday, 2019): the baseline is used as a trial-level regressor instead of being subtracted, and the corrected epochs are stored in the output dataset, ready for further analyses (BASELINE_REGRESSION can also be called directly with condition labels).
- New 'rr' heart signal: provide beat times detected elsewhere ('beat_latencies', in s) for HEP or HRV features.
- HEP: epochs are now -300 to 600 ms by default ('hep_window'; was -300 to 700 ms), and heartbeats followed by the next one before the epoch end + 50 ms are rejected, so no epoch contains the next QRS (with -300 to 700 ms and a 550 ms limit, 28% of the sample-data epochs did). 'hep_window','adaptive' sets the epoch end from the subject's heart rate, for within-subject analyses. The inter-beat-interval rejection now actually removes trials (all events were epoched before), epochs are time-locked to R-peaks only (not to other events in the file), and to the beats kept after RR cleaning.
- HRSP (heartbeat-related spectral perturbations; was 'HEO'): 5-20 Hz, 5-cycle wavelets (Lee et al., 2024), computed on padded epochs and expressed relative to the whole cardiac cycle (the pre-R-peak baseline overlapped the QRS once smeared by the wavelets).
- PPG: HEPs can be time-locked to the heartbeats rather than to the pulses ('ppg_transit': a delay in ms, or an ECG channel to estimate the pulse arrival time from; see ESTIMATE_PAT). On the sample data, the arrival time was 424 ms and the PPG-based HEP correlated with the ECG-based HEP at r = 0.86 after correction (r = -0.17 before).
- EEG band power: theta, alpha, beta and gamma used wrong frequency edges (bin indices instead of Hz); fixed. The 'individualized' band option works again: alpha is bounded by the individual alpha peak (median across channels), with the theta and beta edges moved accordingly (the band limits used are in frequency.bands). IAF threshold (restingIAF) restored to log10. Alpha asymmetry pairs each left electrode with its mirror position (it could pair C3 with T8), and its 'asy_norm' option now divides alpha by each channel's total power. Entropy features are computed at the actual resampled rate.
- HRV: pNN50 is now in % (was a fraction rounded to one decimal), PRSA acceleration/deceleration capacities follow Bauer et al. (2006) (were the mean NN at the anchors), LF/HF is no longer divided by total power with 'hrv_norm'. Band powers: unitless with the default normalized Lomb-Scargle periodogram (were multiplied by 1e6 and labelled ms^2), in ms^2 with 'LombScargle', 'welch' and 'fft'; the mean NN is removed before Welch/FFT (it inflated LF and HF), the 'fft' method works (it crashed), and recordings shorter than 34 s return NaN band powers instead of crashing.
- ECG/PPG signal quality: windows where the two detectors agree on no beat now count as bad (were ignored); the PPG SQI now covers the whole recording (only the last 30 s were kept).
- Brain-heart coherence: the MVAR model is now fitted to the data before computing the coherence measures (the data were passed as model coefficients).
- rm_heart: heart components are no longer removed before the heart-component step; heart channel removed unless 'keep_heart'; the heart-locked EEG amplitude (cardiac field artifact) is reported before and after removal. The 'boost' option was removed: it did not change the components found and re-referenced the data through the rescaled heart channel.
- Command line: 'icamethod' (or 'ica_method'), 'conf_thresh', and get_RR options are now parsed; 'parpool','off' and 'eeg_features'/'hrv_features','off' work; parallel pool uses the 'Processes' profile and no longer changes your saved parallel settings. Runs under 'matlab -batch' no longer block on dialogs or the progress bar.
- GUI: reference and filter-type choices were mapped to the wrong options in some modes; fixed.

v1.5 (5/2/2024) - METHOD 4 (brain-heart coherence) added

v1.4 (4/1/2024) - publication JoVE (methods 1, 2, 3)

## When using BrainBeats, please cite:

Cannard, C., Wahbeh, H., & Delorme, A. (2024). BrainBeats as an Open-Source EEGLAB Plugin to Jointly Analyze EEG and Cardiovascular Signals. Journal of visualized experiments: JoVE, (206).


## BrainBeats was used and cited in:

Abdollahpour, N., & Artan, N. S. (2025). Significant interactions in infant operculum regions when exposed to a bilingual environment: a resting-state fNIRS study. Neurophotonics, 12(4), 045012-045012.

Abdullah, J.et al. (2025). Mathematical Decoding of the Correlation Between Different Organs' Activities: A Review. Fractals, volume 33, issue: 09.

Liu, P., Gao, Y., Ballegaard, M., & Puthusserypady, S. (2025, July). A Novel Framework for Real-Time ECG and Blood Pressure Signal Analysis: Enhancing Accuracy and Interaction through Dynamic Quality Evaluation. In Annual International Conference of the IEEE Engineering in Medicine and Biology Society. IEEE Engineering in Medicine and Biology Society. Annual International Conference (Vol. 2025, pp. 1-7).

Carbone, F., Silva, M., Leemann, B., Hund-Georgiadis, M., & Hediger, K. (2025). Registered Report Stage I: Neurological and physiological effects of animal-assisted treatments for patients in a minimally conscious state: a randomized, controlled cross-over study. Neuroscience.

Balasubramanian, K. et al. (2025). Complexity Measures in Biomedical Signal Analysis: A Clinically-Grounded Survey Across EEG, ECG, Intracranial Pressure, and Photoplethysmogram Modalities. IEEE Access.

Remiszewski, M. (2025). Long-term Aerobic Exercise Enhances Interoception and Reduces Symptoms of Depression and Anxiety in Physically Inactive Young Adults: A Randomized Controlled Trial. Psychology of Sport and Exercise, 102939.

Chowdhury, et al. (2025). Neural Signals, Machine Learning, and the Future of Inner Speech Recognition. Frontiers in Human Neuroscience, 19, 1637174.

Naaz, R., & Ahmad, S. (2025). ECG Data Mining Approach for Detection of Arrhythmia Using Machine Learning. In 2025 3rd International Conference on Device Intelligence, Computing and Communication Technologies (DICCT) (pp. 52-57). IEEE.

Cheng, X., Maess, B., & Schirmer, A. (2025). A Pleasure That Lasts: Convergent neural processes underpin comfort with prolonged gentle stroking. Cortex.

Georgaras, E., & Vourvopoulos, A. (2025). Physiological assessment of brain, cardiovascular, and respiratory changes in multimodal motor imagery brain-computer interface training. Research in Biomedical Engineering and Technology, 12(1), 2471680.

Park, S., Ha, J., & Kim, L. (2025). Improving single-trial detection of error-related potentials by considering the effect of heartbeat-evoked potentials in a motor imagery-based brain-computer interface. Computers in Biology and Medicine, 195, 110563.

Perez, T. M., Drake, E., & Sullivan, S. (2024). Assessing central nervous system and peripheral nervous system functioning in resting and non-resting conditions in a healthy adult population: A feasibility study. Chiropractic Journal of Australia (Online), 51(1), 1-32.

Akuthota, S., Rajkumar, K., & Janapati, R. (2024). Intelligent EEG Artifact Removal in Motor ImageryBCI: Synergizing FCIF, FCFBCSP, and Modified DNN with SNR, PSD, and Spectral Coherence Evaluation. In 2024 International Conference on Circuit, Systems and Communication (ICCSC) IEEE.

Ingolfsson et al. (2024). Brainfusenet: Enhancing wearable seizure detection through eeg-ppg-accelerometer sensor fusion and efficient edge deployment. IEEE Transactions on Biomedical Circuits and Systems.

Fields, C., et al. (2024). Search for entanglement between spatially separated Living systems: Experiment design, results, and lessons learned. Biophysica, 4(2), 168-181.

Cannard, C., Delorme, A., & Wahbeh, H. (2024). Identifying HRV and EEG correlates of well-being using ultra-short, portable, and low-cost measurements. bioRxiv, 2024-02.

Arao, H., Suwazono, S., Kimura, A., Asano, H., & Suzuki, H. (2023). Measuring auditory event‐related potentials at the external ear canal: A demonstrative study using a new electrode and error‐feedback paradigm. European Journal of Neuroscience, 58(11), 4310-4327.

Goodwin, A. J., et al. (2023). The truth Hertz—synchronization of electroencephalogram signals with physiological waveforms recorded in an intensive care unit. Physiological Measurement, 44(8), 085002.
