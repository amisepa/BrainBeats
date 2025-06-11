<!-- <p align="center"> -->
# BrainBeats
<!-- </p> -->

<p align="center" width="100%">
    <img width="30%" src="https://github.com/amisepa/BrainBeats/blob/main/brainbeats_logo2.png">
</p>

The BrainBeats toolbox, implemented as an EEGLAB plugin, allows joint processing and analysis of EEG and cardiovascular signals (ECG and PPG) for brain-heart interplay research. Both the general user interface (GUI) and command line are supported (see tutorial). BrainBeats currently supports: 1) Heartbeat-evoked potentials (HEP) and oscillations (HEO); 2) Extraction of EEG and HRV features; 3) Extraction of heart artifacts from EEG signals; 4) brain-heart coherence.

 
## THREE METHODS AVAILABLE

<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/blob/main/figures/diagram.png">
</p>

1) Process EEG data for heartbeat-evoked potentials (HEP) analysis using ECG or PPG signals. Steps include signal processing of EEG and cardiovascular signals, inserting R-peak markers into the EEG data, segmentation around the R-peaks with optimal window length, time-frequency decomposition.


<p align="center">
    Example of HEP at the subject level, obtained from simultaneous EEG-ECG signals
</p>
<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/blob/main/figures/fig11.png"> 
</p>

<p align="center">
    Example of HEP at the subject level, obtained from simultaneous EEG-PPG signals
</p>
<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/blob/main/figures/fig17.png">
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
    <img width="50%" src="https://github.com/amisepa/BrainBeats/blob/main/figures/fig21.png"> 
</p>

<p align="center">
    Example of EEG features extracted from sample dataset
</p>
<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/blob/main/figures/fig22.png"> 
</p>


3) Remove heart components from EEG signals using ICA and ICLabel.
   
<p align="center">
    Example of extraction of cardiovascular components from EEG signals
</p>
<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/blob/main/figures/fig27.png"> 
</p>


4) Compute brain-heart coherence (beta version, please test and give feedback)
   
<p align="center">
    Example of several brain-heart coherence measures computed with BrainBeats from simultaneous EEG and ECG signals
</p>
<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/blob/main/figures/coherence_allfreqs.png"> 
</p>

<p align="center">
    Scalp topography showing scalp regions coherent with ECG signal for each frequency band
</p>
<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/blob/main/figures/coherence_topo.png"> 
</p>

## Requirements

- MATLAB installed (https://www.mathworks.com/downloads)
- EEGLAB installed (https://github.com/sccn/eeglab)
- Some data containing EEG and cardiovascular signals (ECG or PPG) within the same file (i.e. recorded simultaneously).
  Or use the tutorial dataset provided in this repository located in the "sample_data" folder. Source: sub-32 in https://nemar.org/dataexplorer/detail?dataset_id=ds003838

## Step-by-step tutorial

See our publication for a step-by-step tutorial using the sample dataset: https://www.jove.com/t/65829/brainbeats-as-an-open-source-eeglab-plugin-to-jointly-analyze-eeg

Full-text preprint: https://www.biorxiv.org/content/10.1101/2023.06.01.543272v3.full

## Version history

v1.5 (5/2/2024) - METHOD 4 (brain-heart coherence) added

v1.4 (4/1/2024) - publication JoVE (methods 1, 2, 3)

## Reference (please cite)

Cannard, C., Wahbeh, H., & Delorme, A. (2024). BrainBeats as an Open-Source EEGLAB Plugin to Jointly Analyze EEG and Cardiovascular Signals. Journal of visualized experiments: JoVE, (206).


## Cited by

Ingolfsson et al. (2024). Brainfusenet: Enhancing wearable seizure detection through eeg-ppg-accelerometer sensor fusion and efficient edge deployment. IEEE Transactions on Biomedical Circuits and Systems.

Naaz, R., & Ahmad, S. (2025). ECG Data Mining Approach for Detection of Arrhythmia Using Machine Learning. In 2025 3rd International Conference on Device Intelligence, Computing and Communication Technologies (DICCT) (pp. 52-57). IEEE.

Cheng, X., Maess, B., & Schirmer, A. (2025). A Pleasure That Lasts: Convergent neural processes underpin comfort with prolonged gentle stroking. Cortex.

Georgaras, E., & Vourvopoulos, A. (2025). Physiological assessment of brain, cardiovascular, and respiratory changes in multimodal motor imagery brain-computer interface training. Research in Biomedical Engineering and Technology, 12(1), 2471680.

Akuthota, S., Rajkumar, K., & Janapati, R. (2024). Intelligent EEG Artifact Removal in Motor ImageryBCI: Synergizing FCIF, FCFBCSP, and Modified DNN with SNR, PSD, and Spectral Coherence Evaluation. In 2024 International Conference on Circuit, Systems and Communication (ICCSC) IEEE.

Cannard, C., Delorme, A., & Wahbeh, H. (2024). Identifying HRV and EEG correlates of well-being using ultra-short, portable, and low-cost measurements. bioRxiv, 2024-02.

Arao, H., Suwazono, S., Kimura, A., Asano, H., & Suzuki, H. (2023). Measuring auditory event‐related potentials at the external ear canal: A demonstrative study using a new electrode and error‐feedback paradigm. European Journal of Neuroscience, 58(11), 4310-4327.
