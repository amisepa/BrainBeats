clear; close all; clc
eeglab; close; 
mainDir = fileparts(which('eegplugin_BrainBeats.m')); cd(mainDir);

%% METHOD 1: Process file for HEP analysis

EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(mainDir,'sample_data'));

% EEG = brainbeats_process(EEG);  % GUI mode
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'},'clean_rr','pchip','clean_eeg',true, ...
    'parpool',false,'gpu',false,'vis',true,'save',true); 
% pop_eegplot(EEG,1,1,1)

%% METHOD 2: Extract EEG and HRV features

EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(mainDir,'sample_data'));

% [~, Features] = brainbeats_process(EEG);  % GUI mode
[~, Features, com] = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'}, 'clean_rr','pchip','clean_eeg',true,'norm',true,...
    'eeg_features', {'time' 'frequency'}, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
    'gpu',false,'parpool',true,'save',true,'vis',true);

%% METHOD 3: Remove heart components from EEG signals

EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(mainDir,'sample_data'));

% EEG = brainbeats_process(EEG);  % GUI mode
[EEG,~,com] = brainbeats_process(EEG,'analysis','rm_heart','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'save',true,'vis',true);




