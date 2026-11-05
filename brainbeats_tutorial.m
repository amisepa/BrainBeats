%% Welcome to the BrainBeats tutorial for command line use
%
% REQUIREMENTS:
%   1) MATLAB (requires a license)
%
%   2) EEGLAB
%   Download: https://github.com/sccn/eeglab
%   Unzip (or clone) it on your computer and add it to the MATLAB path:
%   Home panel > Set path > Add folder > select the eeglab folder > Save >
%   Close
%
%   3) The BrainBeats plugin
%   Download: https://github.com/amisepa/BrainBeats
%   Type 'eeglab' in MATLAB's command window to open EEGLAB. Go to File >
%   Manage extensions > type 'brainbeats' in the search bar > select it in
%   the list, and click Install. Or, if you use Git, clone the repo in
%   eeglab > plugins.
%
% Sample dataset used for the tutorial:
% Raw 64-channel EEG, ECG, and PPG data during 3.8 minutes of resting state
% with eyes open. This file corresponds to sub-032_task-rest_eeg.set and
% sub-032_task-rest_ecg.set merged, downsampled to 250 Hz to speed up
% processing.
%
% The original files can be downloaded here:
% https://nemar.org/dataexplorer/detail?dataset_id=ds003838
% These data were recorded with a Brain Products actiCHamp at the
% Ural Federal University.
% Original sample rate = 1000 Hz; power line frequency = 50 Hz; Ground = Fpz;
% Ref = FCz.
% Note: we artificially modified channel 10 (TP9) so that it is detected as
% a bad channel, since there were no bad channels in this dataset. We also
% added electrode artifacts at the beginning of the file and some muscle
% artifacts at 3-6 s on temporal channels, to illustrate artifact removal.
% The script used to prepare this file is "tests" > "prep_sampledata.m".
%
% You can run each section one by one by clicking in the section and
% pressing CTRL (Windows) or CMD (Mac) + ENTER.
%
% We hope this tutorial and BrainBeats are useful to you!
%
% Cedric Cannard & Arnaud Delorme, 2023

%% Open EEGLAB and get the path to the plugin automatically

clear; close all; clc

% Launch EEGLAB (adds its paths and plugins) and close its GUI
eeglab; close;

% Path to the BrainBeats plugin folder (contains the sample_data folder)
main_path = fileparts(which('brainbeats_process.m'));

% Go to the plugin directory
cd(main_path)

%% METHOD 1: Heartbeat-evoked potentials (HEP) and heartbeat-related spectral perturbations (HRSP)

% Load the sample dataset into EEGLAB
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
% Its HEP shows the cardiac field artifact but little after it. For a
% clearer heartbeat-evoked response (the HEP figures of the README), load
% sample_data/sub-047_rest.set instead (another participant of the same
% study, 5.3 min).

% Process file for HEP analysis using default parameters except for:
%   - 'analysis' set to 'hep' (type of analysis)
%   - 'heart_signal' set to 'ECG' (type of heart signal)
%   - 'heart_channels' set to {'ECG'} (name of the ECG electrode)
%   - 'clean_eeg' set to true to preprocess the EEG data with default
%       parameters (0.5-30 Hz filter, average reference, bad channels, bad
%       epochs, ICA + ICLabel; the ICA is fitted on a 1-Hz high-passed copy)
%   - 'ica_method' set to 1 (Picard, fast) instead of 2 (Infomax, default)
%   - 'hep_surrogates' set to 100 to test the heartbeat locking of the HEP,
%       HRSP and HRPC against 100 surrogate heartbeat trains (beats at
%       random cardiac phases)
%   - 'keep_heart' set to true to keep the heart channel in the output
% Epochs span -300 to 600 ms around the R-peaks by default ('hep_window'),
% and heartbeats followed by the next one within 650 ms are rejected so no
% epoch contains the next QRS. For within-subject analyses, 'hep_window',
% 'adaptive' sets the epoch end from the subject's heart rate instead.
% No baseline correction by default; 'hep_baseline','regression' applies a
% regression-based baseline correction (Alday, 2019; window [-150 -50] ms).
% The HEP (time domain), HRSP (time-frequency power) and HRPC (phase:
% pairwise phase consistency across heartbeats), with the surrogate
% control, are computed on all channels (EEG.brainbeats.hrsp, .surrogate),
% and plotted for the average of the frontocentral channels ('hep_roi').
% Other options: 'clean_method','gedai' (GEDAI instead of ICA; the plugin
% is installed if needed) and 'ref','csd' (surface Laplacian).
% Note: the toolbox detects the PPG channel as a non-EEG channel and asks
% to remove it. This is expected: it does not process ECG and PPG at the
% same time.
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true,'ica_method',1, ...
    'hep_surrogates',100,'keep_heart',true);

%% Same, on the independent components (ICs) too
% 'hep_level' sets where the measures are computed: 'channels' (default),
% 'ics' or 'both'. On ICs: the HEP, HRSP, HRPC (and surrogate control) of
% every IC of the cleaning ICA, stored with their scalp maps and ICLabel
% classes in EEG.brainbeats.ics. No IC is selected: heartbeat-locked ICs can
% hold cardiac field artifact, so which ones to analyze is left to you.
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'PPG'});
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true,'ica_method',1, ...
    'hep_level','both','hep_surrogates',100);

%% Same as above but using the PPG signal and adjusting some parameters
%  We change these parameters for demonstration only: default parameters
%  are recommended unless you have a reason to change them.

% Load the file again since it was modified by the previous call
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));

% Here we change the following parameters:
%   - 'heart_signal' set to 'PPG' (signal type)
%   - 'heart_channels' set to {'PPG'} (electrode name)
%   - 'linenoise' set to 50 (power line frequency in Hz)
%   - 'ref' set to 'infinity' to rereference EEG data to infinity instead
%       of common average (default) or 'csd' for current source density
%       transformation (surface Laplacian)
%   - 'highpass' set to 1 to remove EEG frequencies < 1 Hz (default 0.5)
%   - 'lowpass' set to 20 to remove EEG frequencies > 20 Hz
%   - 'filttype' set to 'causal' to use a causal minimum-phase FIR filter
%       instead of the default noncausal zero-phase FIR filter (useful when
%       examining the pre-heartbeat period)
%   - 'detectMethod' set to 'median' to detect and remove bad EEG epochs
%       instead of the default 'grubbs'
%   - 'icamethod' set to 1 (Picard, fast) instead of 2 (Infomax, default)
%       or 3 (replicable Infomax, very slow)
%   - 'save' set to false to not save the final 'filename_HEP.set' file
%   - 'vis_cleaning' set to true to visualize preprocessing plots
%   - 'vis_outputs' set to true to visualize the final outputs
% PPG pulses reach the sensor after the heartbeat (pulse arrival time, PAT:
% 92-112 ms to the pulse onsets on our test data), so the PPG beats are
% shifted back by the PAT ('ppg_transit','auto', default): measured here from
% the ECG channel of the file. Without an ECG, give the delay in ms, or 'eeg'
% to estimate it from the cardiac field artifact of the EEG (104 ms on these
% data, vs 112 ms with the ECG; check its figure: on 6 test subjects it was
% right on 3, wrong on 1, and not clear enough on 2).
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'clean_eeg',true,'linenoise',50, ...
    'ref','infinity','highpass',1,'lowpass',20,'filttype','causal', ...
    'detectMethod','median','icamethod',1, ...
    'save',false,'vis_cleaning',true,'vis_outputs',true);

%% METHOD 2: Extract EEG and HRV features using default parameters

% Load the same raw dataset again
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));

% Launch with default parameters, except:
%   - 'analysis' set to 'features' to extract EEG and HRV features
%   - 'clean_eeg' set to true and 'linenoise' to 50 (Hz)
%   - 'parpool' set to true to speed up the EEG features with parallel computing
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true,'linenoise',50,'parpool',true);

% All features are in EEG.brainbeats.features and, when 'save' is true
% (default), in filename_features.mat next to the .set file loaded in EEGLAB.

% You can replot features using:
% params.chanlocs = EEG.chanlocs;
% plot_features(EEG.brainbeats.features,params)

%% Same but this time, we use the PPG signal and modify some parameters
% Again, this is for illustration only; we recommend default parameters.

% We modify these parameters:
%   - 'clean_eeg' set to false to skip EEG preprocessing
%   - 'hrv_features' set to {'time' 'frequency' 'nonlinear'}
%   - 'eeg_features' set to {'time' 'frequency'} to skip the (slow)
%       nonlinear features
%   - 'hrv_spec' set to 'welch' to use Welch's method on the NN series
%       resampled at 7 Hz instead of the default Lomb-Scargle periodogram
%   - 'eeg_norm' set to 0 to NOT convert PSD to decibels (dB; default = 1)
%   - 'parpool' set to 'on' to use parallel computing
%   - 'save' set to false to not save outputs
%   - 'vis_cleaning' set to false since we already saw them above
%   - 'vis_outputs' set to true to see the outputs
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'ECG'});  % remove ECG channel to avoid warning
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'clean_eeg',false, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'},'hrv_spec','welch', ...
    'eeg_features', {'time' 'frequency'},'eeg_norm',0,...
    'parpool','on','save',false,'vis_cleaning',false,'vis_outputs',true);


%% HRV features only
% To turn OFF all EEG operations, set the input 'eeg' to 'off'.

% ECG
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'PPG'});  % remove PPG channel
EEG = brainbeats_process(EEG,'analysis','features','eeg','off',...
    'heart_signal','ECG', 'heart_channels',{'ECG'}, ...
    'hrv_features',{'time' 'frequency' 'nonlinear'},...
    'vis_cleaning',true,'vis_outputs',true);

% Same but without any EEG data (we remove all EEG channels)
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'channel',{'ECG'}); 
EEG = brainbeats_process(EEG,'analysis','features','eeg','off',...
    'heart_signal','ECG','heart_channels',{'ECG'});

% Same but for PPG (no EEG)
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'channel',{'PPG'}); 
EEG = brainbeats_process(EEG,'analysis','features','eeg','off', ...
    'heart_signal','PPG','heart_channels',{'PPG'});

% Preprocessing outputs can be found in:
EEG.brainbeats.preprocessings

% HRV features can be found in:
EEG.brainbeats.features.HRV

%% Extract EEG features only
% To turn OFF all heart operations, set the input 'heart_signal' to 'off'.

EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'PPG' 'ECG'});  % remove heart channels
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','off', ...
    'eeg_features',{'time' 'frequency','nonlinear'},'clean_eeg',1,'linenoise',50, ...
    'ref','infinity','ica_method',1,'parpool',1,'vis_cleaning',1,'vis_outputs',1,'save',1);

% Preprocessing outputs can be found in:
EEG.brainbeats.preprocessings

% EEG features can be found in:
EEG.brainbeats.features.EEG


%% METHOD 3: Brain-heart coupling (command line only)
% Coupling between the slow fluctuations of EEG band power (per region) and
% heart rate variability, with respiration as a conditioning variable:
% Granger causality in both directions (conditional and pairwise), transfer
% entropy, coherence and partial coherence, directed coherence and GPDC.
% See COMPUTE_COUPLING. Results: EEG.brainbeats.coupling.table.

% ECG, respiration derived from the ECG (no respiration channel here)
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'PPG'});  % remove PPG channel to avoid warning
EEG = brainbeats_process(EEG,'analysis','coupling','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'coupling_resp','ecg','clean_eeg',true,'linenoise',50, ...
    'parpool',false,'vis_outputs',true);
EEG.brainbeats.coupling.table

% PPG (no ECG-derived respiration: HF coupling may be the respiratory drive)
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'ECG'});  % remove ECG channel to avoid warning
EEG = brainbeats_process(EEG,'analysis','coupling','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'coupling_resp','none','clean_eeg',true,'linenoise',50, ...
    'parpool',false,'vis_outputs',true);


%% To launch the BrainBeats window via command line (or EEGLAB menu 'BrainBeats')

EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
[EEG, com] = pop_brainbeats(EEG);

% Display <com> after the run to get the command line equivalent of the
% parameters selected in the GUI
com

