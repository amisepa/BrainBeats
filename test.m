clear; close all; clc
cd('C:\Users\Tracy\Documents\MATLAB\BrainBeats');
% dataDir = 'C:\Users\Tracy\Documents\MATLAB\BrainBeats\sample_data';
eeglab; close;
dataDir = fileparts(which('pop_BrainBeats.m'));

%% MODE 1: Remove heart components from sample_data1

EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(dataDir,'sample_data'));
% EEG = pop_BrainBeats(EEG);  % GUI mode
pop_BrainBeats(EEG,'analysis','rm_heart','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'vis',true);

%% MODE 2: Run HEP on sample_data2

EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(dataDir,'sample_data'));
% EEG = pop_BrainBeats(EEG);  % GUI mode
pop_BrainBeats(EEG,'analysis','hep','heart_signal',{'ECG'}, ...
    'heart_channels',{'EXG5' 'EXG6'},'clean_eeg',true,'vis',true); 

%% MODE 3: Feature-based

EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(dataDir,'sample_data'));
% Features = pop_BrainBeats(EEG);  % GUI mode
Features = pop_BrainBeats(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'EXG5' 'EXG6'}, 'clean_eeg',true, ...
    'eeg_features', {'frequency' 'nonlinear'}, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'},'vis',true);

%% MODE 3: Feature-based (Arno's file)

EEG = pop_loadset('filename','sub-024_task-Post1_eeg.set','filepath','C:\Users\Tracy\Downloads');
Features = pop_BrainBeats(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG'}, 'clean_eeg',false, ...
    'eeg_features', {'frequency' 'nonlinear'}, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'},'vis',true);


%% Save figures (edit name)

exportgraphics(gcf, fullfile('figures','HRV_power_entropy.png'),'Resolution',300)

exportgraphics(gcf, fullfile('figures','GUI_mode3.png'),'Resolution',300)






%% MODE 2: Run HEP anlaysis on group

% Continuous data: mindwandering vs trance
% EEG = pop_biosig('G:\Shared drives\Grants\Post Award Grants\(736) Bial Full-trance 2017\Research\Data\EEG\BDF_files\subj06_1.bdf');
% EEG = pop_select(EEG, 'rmchannel',{'EXG1','EXG2','EXG3','EXG4','EXG7','EXG8','GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp'});
for iSub = 1%1:13
    for iCond = 1%:2
        switch iCond
            case 1
                filename = sprintf('sub-%2.2d_mindwandering.set',iSub);
            case 2
                filename = sprintf('sub-%2.2d_trance.set',iSub);
        end
        EEG = pop_loadset('filename',filename,'filepath',dataDir);
        pop_BrainBeats(EEG,'analysis','hep','heart_signal',{'ECG'}, ...
            'heart_channels',{'EXG5' 'EXG6'},'clean_eeg',true,'vis',true); 

    end
end

% Muse EEG + ECG
% EEG = pop_loadset('filename','0b17e1cada_clean.set','filepath','G:\Shared drives\Science\IDL\5. DATA\muse\eeg\eeg_ecg_clean');

% Muse EEG + PPG
% EEG = import_edf('G:\Shared drives\Science\IDL\5. DATA\muse\eeg\edf_museS\2022-08-01T10_49_02-07_00_6002-PUYU-5DC8_eeg.edf', 1);

% BDF ERP
% EEG = pop_biosig('G:\Shared drives\Grants\Post Award Grants\CLOSED PROJECTS\(737) Bial Mediumship 2017\Research\Data\MD_EEG\MD02\MD02.bdf');
% EEG = pop_select(EEG, 'rmchannel',{'EXG3','EXG4','EXG5','EXG6','EXG7','EXG8','GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp'});


%% Run plugin

% MODE 1, 2, or 3 using GUI

% MODE 2: HEP
pop_BrainBeats(EEG,'heart_signal','ECG','heart_channels',{'EXG5' 'EXG6'}, ...
    'analysis','hep','vis',true);

% MODE 3: Feature-based
outputs = pop_BrainBeats(EEG,'heart_signal','ECG','heart_channels',{'EXG5' 'EXG6'}, ...
    'analysis','features', 'eeg_features', {'frequency' 'nonlinear'}, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'}, 'vis',true);

% MODE 3 (precising which features to compute)
% [outputs, com] = pop_BrainBeats(EEG,'heart_signal','ECG','heart_channels',
% {'EXG1' 'EXG2'},'analysis','continuous','eeg_features',{'PSD' 'IAF' 'ASY' 'COH' 'ENT'}, ...
%     'hrv_features', {'SDNN' 'RMSSD' 'pNN50' 'ULF' 'VLF' 'LF' 'HF' 'LF/HF' 'Total' 'ENT' 'FRAC' 'DFA'},...
%     'vis',true);

