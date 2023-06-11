# BrainBeats

The BrainBeats toolbox, implemented as an EEGLAB plugin, allows joint processing and analysis of EEG and cardiovascular signals (ECG/PPG). Both general user interface (GUI) and command line are supported (see wiki). 

It has 3 main modes: 

  1) Remove heart components from EEG signals using ICA and ICLabel; 

  2) Perform heartbeat-evoked potentials (HEP) analysis, including signal processing, segmentation, time-frequency decomposition, and advanced statistical analysis using hierarchichal linear modeling (LIMO plugin). Statistical analyses can be done using the EEGLAB STUDY mode, using hierarchichal linear modeling and advanced corrections for type 1 error using the LIMO-EEG plugin. 

  3) Extract EEG and HRV features from continuous data in the time, frequency, and nonlinear domains. 
     - HRV time domain: NN statistics, SDNN, RMSSD, pNN50.
     - HRV frequency domain: VLF-power, ULF-power, LF-power, HF-power, LF/HF ratio, Total power. Default method is the normalized Lombscargle-periodogram, but Welch and FFT with resampling are also available. Eahc band-power can be normalized to total power.
     - HRV nonlinear domain: Poincare SD1, SD2, SD1/SD2, fuzzy entropy, multiscale fuzzy entropy, PRSA-AC, PRSA-DC. 

     - EEG time domain: RMS, mode, variance, skewness, kurtosis, interquartile range
     - EEG frequency domain: band power (delta, theta, alpha, gamma), individual alpha frequency (IAF, using center of gravity), alpha asymmetry (all possible pairs, normalized), EEG coherence (all pairs except between enighbors; see Nunez et al. 2001). 
     - EEG nonlinear domain: fuzzy entropy, multiscale fuzzy entropy on all channels

For 2) and 3), the toolbox implements automated cleaning of EEG signals (although pre-processed data is recommended for fine-tuning of parameters), bandpass filter of ECG signal, detection of R peaks, a signal quality index of the RR series (SQI, developed by Vest et al. 2017), and removal/interpolation of RR artifacts to obtain NN intervals (several interpolation methods are available). 

For a step-by-step tutorial using the sample data, please see our preprint: https://www.biorxiv.org/content/10.1101/2023.06.01.543272v1.full.pdf
