clear; close all; clc
eeglab; close; 
mainDir = fileparts(which('eegplugin_BrainBeats.m')); cd(mainDir);

%% METHOD 1: Process file for HEP analysis

EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(mainDir,'sample_data'));
% chans = contains({EEG.chanlocs.labels}, {'PO4'});
% EEG.data(chans,:) = EEG.data(chans,:).*5;
% EEG.event = []; EEG.urevent = []; EEG = eeg_checkset(EEG);

% EEG = brainbeats_process(EEG);  % GUI mode
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'},'clean_rr','pchip','clean_eeg',true, ...
    'parpool',false,'gpu',false,'vis',true,'save',true); 
% pop_eegplot(EEG,1,1,1)

%% METHOD 2: Extract EEG and HRV features

EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(mainDir,'sample_data'));
% [~, Features] = brainbeats_process(EEG);  % GUI mode
[~, Features] = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'}, 'clean_rr','pchip','clean_eeg',true,'norm',true,...
    'eeg_features', {'time' 'frequency' 'nonlinear'}, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
    'gpu',false,'parpool',true,'save',true,'vis',true);

%% METHOD 3: Remove heart components from EEG signals

EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG);  % GUI mode
EEG = brainbeats_process(EEG,'analysis','rm_heart','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'save',true,'vis',true);

%% METHOD 2: GROUP STATS

% Download .zip dataset from: https://nemar.org/dataexplorer/detail?dataset_id=ds003690

clear; close all; clc
eeglab; close; 
chanlocpath = fileparts(which('dipfitdefs.m'));

% Import participants data
dataDir = 'C:\Users\Tracy\Downloads\ds003690';
tmp = readtable(fullfile(dataDir, "participants.tsv"),"FileType","text",'Delimiter', '\t');
ids = tmp.participant_id;
age = tmp.age;
sex = tmp.sex;
grp = tmp.group;
edu = tmp.education_years;
clear tmp

% Process "passive" data and extract features
for iFile = 1:length(ids)
    filepath = fullfile(dataDir,ids{iFile}, 'eeg');
    filename = sprintf('%s_task-passive_run-1_eeg.set',ids{iFile});
    EEG = pop_loadset('filename', filename,'filepath',filepath);    
    % EEG = pop_select(EEG,'chantype', {'ECG','EEG'});
    EEG = pop_select(EEG,'nochannel', {'CB1' 'CB2' 'VEO' 'HEO' 'R-Dia-X-(mm)' 'R-Dia-Y-(mm)'});
    figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    EEG = pop_chanedit(EEG,'rplurchanloc',1,'lookup',fullfile(chanlocpath,'standard_BEM','elec','standard_1005.elc'));
    % pop_eegplot(EEG,1,1,1);
    EEG.ref = 'unknown'; EEG = eeg_checkset(EEG); % force reference to unknown to re-reference to Infinity
    EEG = pop_resample(EEG,256); % downsample to 256 HZ to increase speed
    [~, Features] = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'EKG'},'clean_eeg',true,'norm',true,...
    'eeg_features', {'time' 'frequency'}, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
    'gpu',false,'parpool',true,'save',true,'vis',false);
end
disp('Done! Now, use brainbeats_analyze for group analysis'); gong

%% Muse EEG + PPG

% EEG = import_edf('G:\Shared drives\Science\IDL\5. DATA\muse\eeg\edf_museS\2022-08-01T10_49_02-07_00_6002-PUYU-5DC8_eeg.edf', 1);

%% Save figures for paper (edit name)

exportgraphics(gcf, fullfile('figures','method3_clean-signal.png'),'Resolution',300)
exportgraphics(gcf, fullfile('figures','method3_clean-signal.eps'),'Resolution',300)
% print(gcf,fullfile('figures','hep_bad-channels.png'),'-dpng','-r300');     %300 dpi .png
% print(gcf,fullfile('figures','hep_bad-channels.pdf'),'-dpdf','-r300');     %300 dpi .pdf
% print(gcf,fullfile('figures','hep_bad-channels.epsc'),'-depsc','-r300');   %300 dpi .epsc
