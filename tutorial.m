clear; close all; clc
mainDir = fileparts(which('eegplugin_BrainBeats.m'));
cd(mainDir);
eeglab; close;
 
%% METHOD 1: Process file for HEP analysis

EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(mainDir,'sample_data'));
% chans = contains({EEG.chanlocs.labels}, {'PO4'});
% EEG.data(chans,:) = EEG.data(chans,:).*5;
% EEG.event = []; EEG.urevent = []; EEG = eeg_checkset(EEG);

% EEG = brainbeats_process(EEG);  % GUI mode
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'},'clean_eeg',true, ...
    'save',false,'gpu',false,'vis',true); 
% pop_eegplot(EEG,1,1,1)

%% METHOD 2: Extract EEG and HRV features

EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(mainDir,'sample_data'));
% Features = brainbeats_process(EEG);  % GUI mode
% Features = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
%     'heart_channels',{'EXG5' 'EXG6'}, 'clean_eeg',true, ...
%     'eeg_features', {'frequency' 'nonlinear'}, ...
%     'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
%     'gpu',false,'parpool',true,'vis',true);

[~, Features] = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'}, 'clean_eeg',true, ...
    'eeg_features', {'time'},'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
    'gpu',false,'parpool',true,'vis',true);

%% METHOD 2: GROUP STATS




%% METHOD 3: Remove heart components from EEG signals

% EEG = brainbeats_process(EEG);  % GUI mode
EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','rm_heart','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'vis',true);

%% Muse EEG + PPG

% EEG = import_edf('G:\Shared drives\Science\IDL\5. DATA\muse\eeg\edf_museS\2022-08-01T10_49_02-07_00_6002-PUYU-5DC8_eeg.edf', 1);

%% Save figures for paper (edit name)

exportgraphics(gcf, fullfile('figures','method1_ecg-rr.png'),'Resolution',300)
exportgraphics(gcf, fullfile('figures','method3_ecg-rr.eps'),'Resolution',300)
% print(gcf,fullfile('figures','hep_bad-channels.png'),'-dpng','-r300');     %300 dpi .png
% print(gcf,fullfile('figures','hep_bad-channels.pdf'),'-dpdf','-r300');     %300 dpi .pdf
% print(gcf,fullfile('figures','hep_bad-channels.epsc'),'-depsc','-r300');   %300 dpi .epsc
