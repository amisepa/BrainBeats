# BrainBeats

The BrainBeats toolbox, implemented as an EEGLAB plugin, allows joint processing and analysis of EEG and cardiovascular signals (ECG/PPG). Both general user interface (GUI) and command line are supported (see wiki). 

It has 3 main modes: 

  1) Remove heart components from EEG signals using ICA and ICLabel; 

  2) Perform heartbeat-evoked potentials (HEP) analysis, including signal processing, segmentation, time-frequency decomposition, and advanced statistical analysis using hierarchichal linear modeling (LIMO plugin). 

  3) Extract EEG and HRV features from continuous data, in the time, frequency, and nonlinear domains. 

For 2) and 3), the toolbox implements automated cleaning of EEG signals, detection of R peaks, a signal quality index of the RR series (SQI, developed by Vest et al. 2017), and removal/interpolation of RR artifacts to obtain NN intervals. 
