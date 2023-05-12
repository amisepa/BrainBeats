clear; close all; clc
cd('C:\Users\Tracy\Documents\MATLAB\BrainBeats');
% dataDir = 'C:\Users\Tracy\Documents\MATLAB\BrainBeats\sample_data';
eeglab; close;
dataDir = fileparts(which('pop_BrainBeats.m'));

% Muse EEG + ECG
% EEG = pop_loadset('filename','0b17e1cada_clean.set','filepath','G:\Shared drives\Science\IDL\5. DATA\muse\eeg\eeg_ecg_clean');

% Muse EEG + PPG
% EEG = import_edf('G:\Shared drives\Science\IDL\5. DATA\muse\eeg\edf_museS\2022-08-01T10_49_02-07_00_6002-PUYU-5DC8_eeg.edf', 1);

% BDF ERP
% EEG = pop_biosig('G:\Shared drives\Grants\Post Award Grants\CLOSED PROJECTS\(737) Bial Mediumship 2017\Research\Data\MD_EEG\MD02\MD02.bdf');
% EEG = pop_select(EEG, 'rmchannel',{'EXG3','EXG4','EXG5','EXG6','EXG7','EXG8','GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp'});

%% MODE 1: Remove heart components from sample_data1

EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(dataDir,'sample_data'));
% EEG = pop_BrainBeats(EEG);  % GUI mode
pop_BrainBeats(EEG,'analysis','rm_heart','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'vis',true);

%% MODE 2: Run HEP on sample_data2

EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(dataDir,'sample_data'));
% EEG = pop_BrainBeats(EEG);  % GUI mode
EEG = pop_BrainBeats(EEG,'analysis','hep','heart_signal',{'ECG'}, ...
    'heart_channels',{'EXG5' 'EXG6'},'clean_eeg',true,'gpu',true,'vis',true); 

%% MODE 3: Feature-based

EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(dataDir,'sample_data'));
% EEG.data(30,:) = EEG.data(30,end:-1:1).*10;
% EEG.data(50,:) = EEG.data(50,end:-1:1).*10;
% Features = pop_BrainBeats(EEG);  % GUI mode
[EEG, Features] = pop_BrainBeats(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'EXG5' 'EXG6'}, 'clean_eeg',true, ...
    'eeg_features', {'frequency' 'nonlinear'}, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
    'gpu',true,'parpool',true,'vis',true);

%% MODE 3: Feature-based (Arno's file)

% EEG = pop_loadset('filename','sub-024_task-Post1_eeg.set','filepath','C:\Users\Tracy\Downloads');
% [EEG, Features] = pop_BrainBeats(EEG,'analysis','features','heart_signal','ECG', ...
%     'heart_channels',{'ECG'}, 'clean_eeg',true, ...
%     'eeg_features', {'frequency' 'nonlinear'}, ...
%     'hrv_features', {'time' 'frequency' 'nonlinear'},'vis',true);


%% Save figures (edit name)

exportgraphics(gcf, fullfile('figures','HRV_power_entropy.png'),'Resolution',300)



%% Run HEP analysis on group: Mindwandering vs Trance

dataDir = 'C:\Users\Tracy\Desktop\trance_sample_data';

commands = {};
count = 0;
for iSub = 1%1:13
    for iCond = 1%:2
        if iCond == 1
            filename = sprintf('sub-%2.2d_mindwandering.set',iSub);
        else
            filename = sprintf('sub-%2.2d_trance.set',iSub);
        end
        EEG = pop_loadset('filename',filename,'filepath',dataDir);

        EEG = pop_BrainBeats(EEG,'analysis','hep','heart_signal',{'ECG'}, ...
            'heart_channels',{'EXG5' 'EXG6'},'clean_eeg',true,'gpu',true,'vis',true); 
        EEG.subject = sprintf('sub-%2.2d',iSub);
        EEG.session = 1;
        if iCond == 1
            EEG.condition = {'mindwandering'};
        else
            EEG.condition = {'trance'};
        end

        EEG.saved = 'no';
        pop_saveset(EEG, 'filepath', EEG.filepath, 'filename', filename);

        count = count+1;
        commands = {commands{:} 'index' count 'load' fullfile(EEG.filepath, filename)};

    end
end

% Create study
[STUDY, ALLEEG] = std_editset([],[],'name','Trance','task','Trance',...
        'filename','trance.study','resave', 'on', ...
        'filepath', dataDir,'commands',commands);
EEG = ALLEEG; CURRENTSET = 1:length(ALLEEG); CURRENTSTUDY = 1;

% % plot results
% STUDY = pop_specparams(STUDY, 'topofreq',[9 11] );
% STUDY = pop_statparams(STUDY, 'condstats','on','mcorrect','fdr');
% STUDY = std_specplot(STUDY,ALLEEG,'channels',{'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz','Fpz','Fp2','AF8','AF4','Afz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'}, 'design', 1);
% print -djpeg channeling_10Hz.jpg
% STUDY = pop_specparams(STUDY, 'topofreq',[4 6] );
% STUDY = std_specplot(STUDY,ALLEEG,'channels',{'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz','Fpz','Fp2','AF8','AF4','Afz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'}, 'design', 1);
% print -djpeg channeling_5Hz.jpg
% STUDY = pop_specparams(STUDY, 'topofreq',[18 22] );
% STUDY = std_specplot(STUDY,ALLEEG,'channels',{'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz','Fpz','Fp2','AF8','AF4','Afz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'}, 'design', 1);
% print -djpeg channeling_20Hz.jpg
% 
% gong