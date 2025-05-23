%% Welcome to the BrainBeats tutorial for command line use
% 
% REQUIREMENTS:
%   1) MATLAB (requires a license) or Octave installed
% 
%   2) Install EEGLAB
%   Download here: https://github.com/sccn/eeglab
%   Unzip (or clone) the file on your computer and add the path to MATLAB: 
%   Home panel > Set path > Add folder > select the eeglab folder > Save >
%   Close
% 
%   3) Install the BrainBeats plugin
%   Download here: https://github.com/amisepa/BrainBeats
%   Type 'eeglab' in MATLAB's command window to open EEGLAB. Go to File >
%   Manage extensions > type 'brainbeats' in the search bar > select in the
%   list area, and click Install. Or, if you use Git, simply clone the repo 
%   in eeglab > plugins on your computer. 
%  
% Sample dataset used for the tutorial:
% Raw 64-channel EEG, ECG, and PPG data during 3.8 minutes of resting state
% with eyes opened. This file corresponds to sub-032_task-rest_eeg.set and
% sub-032_task-rest_ecg.set merged, downsampled to 250 hz to accelerate 
% operations. 
% 
% The original files can be downloaded here:
% https://nemar.org/dataexplorer/detail?dataset_id=ds003838
% These data were recorded with a Brain Products actiCHamp at the 
% Ural Federal University. 
% Original sample rate = 1000 Hz; power line frequency = 50 Hz; Ground = Fpz; 
% Ref = FCz. 
% Note: we artificially modified channel 10 (TP9) to be detected as a bad 
% channel for demonstration purposes since there were no bad channels in this
% dataset. And we artifically added electrode artifacts in the beginning of
% the file, and some muscle artifacts at 3-6 s on temporal channels, for
% illustration of artifact removal. The script used to prepare this file can 
% be found in "functions" > "prep_sampledata.m"
% 
% You can launch each section one by one by clicking in the section and 
% pressing CTRL (for Windows) or CMD (Mac) + ENTER
% 
% We hope you this tutorial and BrainBeats are useful to you too! 
% 
% Cedric Cannard & Arnaud Delorme, 2023

%% Open EEGLAB and get the path to the plugin automatically

clear; close all; clc

% Launch EEGLAB without the GUI
eeglab; close;

% Get the path to the EEGLAB plugin which contains a folder with the sample data
main_path = fileparts(which('brainbeats_process.m'));

% Go to the plugin directory
cd(main_path)

%% METHOD 1: Heartbeat-evoked potentials (HEP) and oscillations (HEO)

% Load the sample dataset into EEGLAB
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));

% Process file for HEP analysis using default parameters except for:
%   - selecting the type of analysis: 'hep'
%   - selecting the type of heart signal: 'ECG'
%   - selecting the name of the ECG electrodes: 'ECG'
%   - 'clean_eeg' set to to 'true' preprocess the EEG data with default parameters.
% Note: the toolbox automatically detects the undesired PPG channel, which 
% is expected since the toolbox is not designed to run both ECG and PPG at 
% the time.
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true,'ica_method',1,'keep_heart',true);

%% Same as above but using the PPG signal and adjusting some parameters 
%  Note that we are changing these parameters for demonstration only, but
%  default parameters are recommended. These parameters should only be
%  changed if you know why. 

% We need to load the file again since it was modified by our previous call above. 
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));

% Here we change the following parameters:
%   - 'heart_signal' set to 'PPG' (signal type)
%   - 'heart_channels' set to 'PPG' (electrode name)
%   - 'clean_rr' set to 'spline' to interpolate the RR artifacts instead of
%       'pchip' (default)
%   - 'ref' set to 'infinity' to rereference EEG data to infinity instead
%       of average (default)
%   - 'highpass' filter set to .5 to remove EEG frequencies <0.5 hz 
%   - 'lowpass' set to 20 to remove EEG frequencies >20 hz
%   - 'filttype' set to 'causal' to use causal minimum-phase FIR filter
%       instead of the default noncausal zero-phase FIR filter (useful if
%       examining the pre-heartbeat period)
%   - 'detectMethod' set to 'median' to detect and remove bad 
%       EEG epochs instead of the default 'grubbs'. 
%   - 'icamethod' to 1 (fast picard) instead of 2 (Infomax) or 3 (modified
%       Infomax, very slow but replicable)
%   - 'keep_heart' to true to preserve the heart channel in final output
%       (e.g. for visual check of final HEP output).
%   - 'save' set to false to not save the final 'filename_HEP.set' file
%   - 'vis_cleaning' set to true to visualize preprocessing plots
%   - 'vis_outputs' set to true to visualize the final outputs 
% Note: the toolbox automatically detects the undesired ECG channel, 
% which is expected since the toolbox is not designed to run both ECG and PPG at the time.
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'clean_rr','spline','clean_eeg',true,'linenoise',50, ...
    'ref','infinity','highpass',.5,'lowpass',20,'filttype','causal', ...
    'detectMethod','median','icamethod',1, ...
    'save',false,'vis_cleaning',true,'vis_outputs',true);

%% METHOD 2: Extract EEG and HRV features using default parameters

% Load the same raw dataset again
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));

% Launch with default parameters
% Note that 'analysis' is set to 'features' to extract EEG and HRV features
% parpool set to ON to accelerate computation of EEG features
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true,'linenoise',50,'parpool',true);

% All features can be found in EEG.brainbeats.features or in a .mat file
% saved in the same place as the .set file loaded in EEGLAB
% (filename_features.mat) if "save" input is set to true.

% You can replot features using:
% params.chanlocs = EEG.chanlocs;
% plot_features(EEG.brainbeats.features,params)

%% Same but this time, we use PPG signal and modify some parameters 
% Again, this is for illustration purpose only, we recommend using default
% parameters. 

% We modify these parameters:
%   - 'hrv_features' to {'time' 'frequency' 'nonlinear'} 
%   - 'eeg_features' to {'frequency'} to only compute frequency-domain
%       features
%   - 'hrv_spec' to 'LombScargle' to use the standard version of the
%       algorithm
%   - 'eeg_norm' set to 0 to NOT convert PSD to decibles (db; default = 1) 
%   - 'parpool' set to false to prevent using parallel computing
%   - 'save' to false to prevent saving
%   - 'vis_cleaning' to 'false' since we already saw them above
%   - 'vis_outputs' to 'true' to see the outputs
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'ECG'});  % remove ECG channel to avoid warning
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'clean_eeg',false,'linenoise',50, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'},'hrv_spec','LombScargle', ...
    'eeg_features', {'time' 'frequency'},'eeg_norm',0,...
    'parpool','on','save',false,'vis_cleaning',false,'vis_outputs',true);


%% HRV features only
% To turn OFF all EEG operations, the input 'eeg_features' is set to 'false'.

% ECG
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'PPG'}); 
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
% To turn OFF all HEART operations, the input 'heart_signal' is set to 'off'.

EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'PPG' 'ECG'});  % remove heart channels
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','off', ...
    'eeg_features',{'time' 'frequency','nonlinear'},'clean_eeg',1,'linenoise',50, ...
    'ref','infinity','ica_method',1,'parpool',1,'vis_cleaning',1,'vis_outputs',1,'save',1);

% Preprocessing outputs can be found in:
EEG.brainbeats.preprocessings

% EEG features can be found in:
EEG.brainbeats.features.EEG


%% METHOD 3: Remove heart components from EEG signals
% to avoid preprocessing the whole file again, remove PPG channel and the
% first 10 s that contain simulated artifacts (that were added for illustration). 

EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','rm_heart','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'linenoise',50,'vis_cleaning',false,...
    'conf_thresh',.8,'boost',true,'ica_method',1);  % 'ica_method',1 for fast PICARD alg


%% METHOD 4: Brain-heart coherence (NEW: BETA; command line only)
% This method has not been tested much yet. Please use with caution and
% report any errors at: https://github.com/amisepa/BrainBeats/issues

% ECG
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'PPG'});  % remove ECG channel to avoid warning
EEG = brainbeats_process(EEG,'analysis','coherence','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',0,'linenoise',50,'ref','infinity','ica_method',1,...
    'parpool',0,'vis_outputs',1);

% PPG
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'ECG'});  % remove ECG channel to avoid warning
EEG = brainbeats_process(EEG,'analysis','coherence','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'coh_signal','hrv','clean_eeg',true,'linenoise',50, ...
    'ref','infinity','ica_method',1,'parpool',false,'vis_outputs',true);


%% To launch the main GUI via command line

EEG = pop_loadset('filename','dataset.set','filepath',fullfile(main_path,'sample_data'));
[EEG, com] = brainbeats_process(EEG);

% Type <com> at the end of the operations to output the command line with
% the parameters that were selected manually in the GUI
com

