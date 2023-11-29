% BrainBeats tutorial
% Launch each section one by one by clicking in the section and pressing: 
% CTRL/CMD + ENTER

clear; close all; clc
eeglab; close; % launch EEGLAB wihtout the GUI

% Find path to the plugin to load sample data files by command line 
addpath('C:\Users\Tracy\Documents\MATLAB\BrainBeats')               %%%%%%%%%%%%%%%%%% TO REMOVE  %%%%%%%%%%%%%%%%%
mainDir = fileparts(which('eegplugin_BrainBeats.m')); 
cd(mainDir); % go to the plugin directory

%% METHOD 1: Heartbeat-evoked potentials (HEP) or oscillations (HEO) 

% Load sample raw file into EEGLAB (64-channel EEG data during 5 minutes of
% resting state, active mindwandering, not preprocessed)
EEG = pop_loadset('filename','data_raw.set','filepath',fullfile(mainDir,'sample_data'));

% Process file for HEP analysis using default parameters except for
% preprocessing the EEG data and visualizing the cleaning performed by the 
% toolbox by setting the inputs 'clean_eeg' and 'vis_cleaning' to 'true'. 
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'},'clean_eeg',true,'vis_cleaning',true); 
% pop_eegplot(EEG,1,1,1);  % to visualize the final output

% Same but adjusting some paramters this time
%   - 'clean_rr' to choose the method to interpolate RR artifacts (default = 'pchip')
%   - 'clean_eeg' to clean the EEG data (true) or not (false) 
%   - 'parpool' to turn ON (true) or OFF (false) the parrallel toolbox
%   - 'gpu' to turn ON (true) or OFF (false) GPU computing
%   - 'save' to save the resulting file (true) or not (false)
%   - 'vis' to visualize outputs with plots (true) or not (false)
% Note that we need to load the file again since it was modified by our 
% previous call above. Here we use spline interpolation of RR artifacts
% (just as an example), no parpool, GPU-computing, no output is saved,
% with visualization  turned ON. 
EEG = pop_loadset('filename','data_raw.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'},'clean_rr','spline','clean_eeg',true, ...
    'parpool',false,'gpu',true,'save',false,'vis_cleaning',true); 


%% METHOD 2: Extract EEG and HRV features

% Load raw dataset again
EEG = pop_loadset('filename','data_raw.set','filepath',fullfile(mainDir,'sample_data'));

% Process file (slightly different approach than for HEP analysis that is
% trial-based, here artifacts are removed with ASR on the continous data). 
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'},'clean_eeg',true,'norm',true,...
    'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
    'eeg_features', {'time' 'frequency' 'nonlinear'}, ...
    'gpu',false,'parpool',true,'save',false,...
    'vis_cleaning',false,'vis_outputs',true);

%% METHOD 3: Remove heart components from EEG signals

EEG = pop_loadset('filename','data_raw.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','rm_heart','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'save',false,'vis',true);

%% To launch the GUI via command line

EEG = brainbeats_process(EEG);



