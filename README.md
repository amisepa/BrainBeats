<!-- <p align="center"> -->
# BrainBeats (Beta)
<!-- </p> -->

<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/blob/v1.3/brainbeats_logo.png">
</p>

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


## Tutorial

See the JoVE preprint for a step-by-step tutorial using the sample dataset: https://www.biorxiv.org/content/10.1101/2023.06.01.543272v1.full.pdf
