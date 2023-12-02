% BrainBeats tutorial
% Launch each section one by one by clicking in the section and pressing:
% CTRL/CMD + ENTER

clear; close all; clc
eeglab; close; % launch EEGLAB wihtout the GUI

% Get the path to the EEGLAB plugin which contains a folder with the sample data
mainDir = fileparts(which('eegplugin_BrainBeats.m'));

% Go to the plugin directory
cd(mainDir); 

%% METHOD 1: Heartbeat-evoked potentials (HEP) and oscillations (HEO)

% Load the dataset1 into EEGLAB. This file contains raw 64-channel EEG,
% ECG, and PPG data during 3.8 minutes of resting state eyes opened. This
% file corresponds to sub-032_task-rest_eeg.set and
% sub-032_task-rest_ecg.set merged, available here:
% https://nemar.org/dataexplorer/detail?dataset_id=ds003838
% note: sub-032_task-rest_ecg.set contains both ECg and PPG channels
EEG = pop_loadset('filename','dataset1.set','filepath',fullfile(mainDir,'sample_data'));
EEG.ref = [];

% Process file for HEP analysis using default parameters except for:
%   - selecting the type of analysis: 'hep'
%   - selecting the type of heart signal: 'ECG'
%   - selecting the name of the ECG electrodes: 'ECG1' and 'ECG2'
%   - preprocessing the raw EEG data by setting the 'clean_eeg' input to 'true'.
% Note that the toolbox detects the PPG channel as undesired automatically, 
% which is to be expected since we cannot run both ECG and PPG at the time.
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true);
% pop_eegplot(EEG,1,1,1);  % to visualize the final output

% Same with PPG
EEG = pop_loadset('filename','dataset1.set','filepath',fullfile(mainDir,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'ECG'}); 
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'eeg',false,'hrv_features',{'time' 'frequency'},...
    'vis_cleaning',true,'vis_outputs',false);

% Same but adjusting some parameters this time
%   - 'clean_rr' set to 'spline' to interpolate the RR artifacts (default = 'pchip')
%   - 'clean_eeg' to preprocess the EEG data (filtering, removal and 
%       interpolation of bad channels, removal of large artifacts, segmentaiton, 
%       removal of bad epochs, removal of artifactual independent
%       components with ICA and ICALabel. 
%   - 'parpool' set to 'false' to turn off parrallel computing (not
%       beneficial here)
%   - 'gpu' set to 'true' to turn off GPU computing (not beneficial here)
%   - 'save' set to 'true' to save the final 'filename_HEP.set' file
%   - 'vis_cleaning' set to 'false' as we already generated all the preprocessing
%       plots above.
%   - 'vis_outputs' set to 'true' to visualize the final output again. 
% Note that we need to load the file again since it was modified by our
% previous call above. 
EEG = pop_loadset('filename','dataset1.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'},'clean_rr','spline','clean_eeg',false, ...
    'parpool',false,'gpu',false,'save',false,'vis_cleaning',false,'vis_outputs',true);




%% new sample data

EEG = pop_loadset('filename','dataset1.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true, ...
    'parpool',false,'gpu',false,'save',true,'vis_cleaning',true,'vis_outputs',true);

EEG = pop_loadset('filename','dataset3_rest.set','filepath',fullfile(mainDir,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'PPG'}); 
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'norm',true,...
    'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
    'eeg_features', {'time' 'frequency' 'nonlinear'}, ...
    'gpu',false,'parpool',true, 'vis_cleaning',true,'vis_outputs',true,...
    'save_plots',true,'save_outputs');

%% METHOD 2: Extract EEG and HRV features

% Load the same raw dataset again
EEG = pop_loadset('filename','dataset1.set','filepath',fullfile(mainDir,'sample_data'));

% This time, note these different inputs:
%   - 'analysis' set to 'features' to extract EEG and HRV features
%   - 'norm' to normalize the frequency domain outputs (type
%       'help get_hrv_features' and 'help get_eeg_features' for more detail)
%   - 'hrv_features' and 'eeg_features' to {'time' 'frequency' 'nonlinear'} 
%       to extract features in all three domains. (type 
%       'help get_hrv_features' and 'help get_eeg_features' for more detail)
%   - 'parpool' to 'true' to use parallel computing.
%   - 'vis_cleaning' to 'true' to see preprocessing steps (slightly different
%       than for HEP). 
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'},'clean_eeg',true,'norm',true,...
    'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
    'eeg_features', {'time' 'frequency' 'nonlinear'}, ...
    'gpu',false,'parpool',true,'save',false,...
    'vis_cleaning',false,'vis_outputs',true);

%% METHOD 3: Remove heart components from EEG signals

EEG = pop_loadset('filename','dataset1_clean.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','rm_heart','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'save',false,'vis',true);

%% To launch the GUI via command line


EEG = brainbeats_process(EEG);




%% Group analysis: HEP



%% Group analysis: Features




