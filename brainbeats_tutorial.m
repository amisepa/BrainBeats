% BrainBeats tutorial
% Launch each section one by one by clicking in the section and pressing: 
% CTRL/CMD + ENTER

clear; close all; clc
eeglab; close; % launch EEGLAB wihtout the GUI

% Find path to the plugin to load sample data files by command line 
mainDir = fileparts(which('eegplugin_BrainBeats.m')); 
cd(mainDir); % go to the plugin directory


%% METHOD 1: Heartbeat-evoked potentials (HEP) or oscillations (HEO) 

% Load sample raw data file into EEGLAB
EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(mainDir,'sample_data'));

% Process file for HEP analysis using default parameters (WARNINGL EEG data 
% must be preprocessed because they will not be preprocessed by default!!!)
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'}); 
% pop_eegplot(EEG,1,1,1);  % to visualize


% Same but adjusting some paramters
%   - 'clean_rr' to choose the method to interpolate RR artifacts (default = 'pchip')
%   - 'clean_eeg' to clean the EEG data (true) or not (false) 
%   - 'parpool' to turn ON (true) or OFF (false) the parrallel toolbox
%   - 'gpu' to turn ON (true) or OFF (false) GPU computing
%   - 'save' to save the resulting file (true) or not (false)
%   - 'vis' to visualize outputs with plots (true) or not (false)
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'},'clean_rr','spline','clean_eeg',false, ...
    'parpool',false,'gpu',true,'save',false,'vis',true); 


%% METHOD 2: Extract EEG and HRV features

EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'}, 'clean_rr','pchip','clean_eeg',true,'norm',true,...
    'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
    'eeg_features', {'time' 'frequency'}, ...
    'gpu',false,'parpool',true,'save',false,'vis',true);

%% METHOD 3: Remove heart components from EEG signals

EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','rm_heart','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'save',false,'vis',true);

%% To launch GUI via command line

EEG = brainbeats_process(EEG);

