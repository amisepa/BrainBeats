<!-- <p align="center"> -->
# BrainBeats
<!-- </p> -->

<p align="center" width="100%">
  <img width="30%" alt="BrainBeats logo"
       src="https://raw.githubusercontent.com/amisepa/BrainBeats/main/brainbeats_logo2.png">
</p>


The BrainBeats toolbox, implemented as an EEGLAB plugin, allows joint processing and analysis of EEG and cardiovascular signals (ECG and PPG) for brain-heart interplay research. Both the general user interface (GUI) and command line are supported (see tutorial). BrainBeats currently supports: 1) Heartbeat-evoked potentials (HEP) and oscillations (HEO); 2) Extraction of EEG and HRV features; 3) Extraction of heart artifacts from EEG signals; 4) brain-heart coherence.

 
## THREE METHODS AVAILABLE

<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/main/figures/diagram.png">
</p>

1) Process EEG data for heartbeat-evoked potentials (HEP) analysis using ECG or PPG signals. Steps include signal processing of EEG and cardiovascular signals, inserting R-peak markers into the EEG data, segmentation around the R-peaks with optimal window length, time-frequency decomposition.


<p align="center">
    Example of HEP at the subject level, obtained from simultaneous EEG-ECG signals
</p>
<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/main/figures/fig11.png"> 
</p>

<p align="center">
    Example of HEP at the subject level, obtained from simultaneous EEG-PPG signals
</p>
<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/main/figures/fig17.png">
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
    <img width="50%" src="https://github.com/amisepa/BrainBeats/main/figures/fig21.png"> 
</p>

<p align="center">
    Example of EEG features extracted from sample dataset
</p>
<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/main/figures/fig22.png"> 
</p>


3) Remove heart components from EEG signals using ICA and ICLabel.
   
<p align="center">
    Example of extraction of cardiovascular components from EEG signals
</p>
<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/main/figures/fig27.png"> 
</p>


4) Compute brain-heart coherence (beta version, please test and give feedback)
   
<p align="center">
    Example of several brain-heart coherence measures computed with BrainBeats from simultaneous EEG and ECG signals
</p>
<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/main/figures/coherence_allfreqs.png"> 
</p>

<p align="center">
    Scalp topography showing scalp regions coherent with ECG signal for each frequency band
</p>
<p align="center" width="100%">
    <img width="50%" src="https://github.com/amisepa/BrainBeats/main/figures/coherence_topo.png"> 
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

## When using BrainBeats, please cite:

Cannard, C., Wahbeh, H., & Delorme, A. (2024). BrainBeats as an Open-Source EEGLAB Plugin to Jointly Analyze EEG and Cardiovascular Signals. Journal of visualized experiments: JoVE, (206).


## BrainBeats was used and cited in:

Carbone, F., Silva, M., Leemann, B., Hund-Georgiadis, M., & Hediger, K. (2025). Registered Report Stage I: Neurological and physiological effects of animal-assisted treatments for patients in a minimally conscious state: a randomized, controlled cross-over study. Neuroscience.

Balasubramanian, K. et al. (2025). Complexity Measures in Biomedical Signal Analysis: A Clinically-Grounded Survey Across EEG, ECG, Intracranial Pressure, and Photoplethysmogram Modalities. IEEE Access.

Park, S. et al. (2025). Improving single-trial detection of error-related potentials by considering the effect of heartbeat-evoked potentials in a motor imagery-based brain-computer interface. Computers in Biology and Medicine, 195, 110563.

Abdullah, J. et al. (2025). Mathematical Decoding of the Correlation Between Different Organs' Activities: A Review. Fractals.

Carbone, F. et al. (2025). Registered Report Stage I: Neurological and physiological effects of animal-assisted treatments for patients in a minimally conscious state: a randomized, controlled cross-over study. Neuroscience.

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
