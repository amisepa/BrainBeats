# BrainBeats


<img src="https://github.com/amisepa/BrainBeats/blob/v1.3/brainbeats_logo.png" width="400">


The BrainBeats toolbox, implemented as an EEGLAB plugin, allows joint processing and analysis of EEG and cardiovascular signals (ECG/PPG). Both the general user interface (GUI) and command line are supported (see tutorial). 

## 3 methods available 

  1) Prepares data for heartbeat-evoked potentials (HEP) analysis, including signal processing, adding R-peaks events into the EEG data, segmentation around the R-peaks, time-frequency decomposition, exporting .set files for statistical analysis. Statistical analyses can then be done using the EEGLAB STUDY mode or hierarchical linear modeling and advanced corrections for the type 1 error with the LIMO-EEG plugin (see LIMO-EEG page for tutorial). 

  2) Extract EEG and HRV features from continuous data in the time, frequency, and nonlinear domains. 
     - HRV time domain: NN statistics, SDNN, RMSSD, pNN50.
     - HRV frequency domain: VLF-power, ULF-power, LF-power, HF-power, LF/HF ratio, Total power. 
     - HRV nonlinear domain: Poincare SD1, SD2, SD1/SD2, fuzzy entropy, multiscale fuzzy entropy, PRSA AC/DC. 

     - EEG time domain: RMS, mode, variance, skewness, kurtosis, interquartile range
     - EEG frequency domain: band power (delta, theta, alpha, gamma), individual alpha frequency (IAF, using the alpha center of gravity), alpha asymmetry (all possible pairs, normalized), EEG coherence.
     - EEG nonlinear domain: fuzzy entropy, multiscale fuzzy entropy on all channels

  3) Remove heart components from EEG signals using ICA and ICLabel; 

## Notes
For HRV frequency metrics, the default method is the normalized Lombscargle-periodogram, as it doesn't require resampling, thus preserving original information. The Welch and FFT methods are available with resampling if desired (option only available via command line currently). Normalization is applied (default) to better deal with irregularly sampled data (common in HRV analysis) during the periodogram estimation by scaling the total power with variance in the time series (this also makes results more comparable across sessions or subjects) and in a 2nd step by dividing each band power by the total power, to provide a more intuitive measure of the relative contribution of each frequency component to the overall power.

For 1) and 2), the toolbox applies a bandpass filter on the ECG signal, detects R peaks, and interpolates (or removes if the user selects this option) RR artifacts to obtain the NN series. For EEG, while preprocessed data, it is always recommended to fine-tune the processing of your data. However, if users choose to use the toolbox to process their EEG data, abnormal channels are detected, removed, and interpolated. For HEP, very bad trials are removed before applying ICA. For features analysis, large artifacts are removed with ASR before performing ICA. The PICARD algorithm is used for fast source separation, using PCA-dimension reduction to account for effective rank deficiency (if present). 


## Tutorial

For a step-by-step tutorial using the sample data, please see our preprint: https://www.biorxiv.org/content/10.1101/2023.06.01.543272v1.full.pdf
