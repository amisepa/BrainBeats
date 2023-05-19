clear; close all; clc
cd('C:\Users\Tracy\Documents\MATLAB\BrainBeats');
% dataDir = 'C:\Users\Tracy\Documents\MATLAB\BrainBeats\sample_data';
eeglab; close;
dataDir = fileparts(which('eegplugin_BrainBeats.m'));

% Muse EEG + ECG
% EEG = pop_loadset('filename','0b17e1cada_clean.set','filepath','G:\Shared drives\Science\IDL\5. DATA\muse\eeg\eeg_ecg_clean');

% Muse EEG + PPG
% EEG = import_edf('G:\Shared drives\Science\IDL\5. DATA\muse\eeg\edf_museS\2022-08-01T10_49_02-07_00_6002-PUYU-5DC8_eeg.edf', 1);

%% MODE 1: Remove heart components from sample_data1

EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(dataDir,'sample_data'));
% EEG = brainbeats_process(EEG);  % GUI mode
EEG = brainbeats_process(EEG,'analysis','rm_heart','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'vis',true);

%% MODE 2: Run HEP on sample_data2

EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(dataDir,'sample_data'));
% EEG = brainbeats_process(EEG);  % GUI mode
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal',{'ECG'}, ...
    'heart_channels',{'ECG1' 'ECG2'},'clean_eeg',true,'gpu',true,'vis',true); 

%% MODE 3: Feature-based

EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(dataDir,'sample_data'));
% EEG.data(30,:) = EEG.data(30,end:-1:1).*10;
% EEG.data(50,:) = EEG.data(50,end:-1:1).*10;
% Features = brainbeats_process(EEG);  % GUI mode
[EEG, Features] = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'EXG5' 'EXG6'}, 'clean_eeg',true, ...
    'eeg_features', {'frequency' 'nonlinear'}, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
    'gpu',false,'parpool',true,'vis',true);

%% MODE 3: Feature-based (Arno's file)

% EEG = pop_loadset('filename','sub-024_task-Post1_eeg.set','filepath','C:\Users\Tracy\Downloads');
% [EEG, Features] = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
%     'heart_channels',{'ECG'}, 'clean_eeg',true, ...
%     'eeg_features', {'frequency' 'nonlinear'}, ...
%     'hrv_features', {'time' 'frequency' 'nonlinear'},'vis',true);


%% Save figures for paper (edit name)

exportgraphics(gcf, fullfile('figures','features_plot1.png'),'Resolution',300)
