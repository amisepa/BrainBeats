%% Prep sample dataset for tutorial

% Load raw data
EEG = pop_loadset('filename','sub-032_task-rest_eeg.set','filepath','C:\\Users\\Tracy\\Documents\\MATLAB\\BrainBeats\\sample_data\\');
CARDIO = pop_loadset('filename','sub-032_task-rest_ecg.set','filepath','C:\\Users\\Tracy\\Documents\\MATLAB\\BrainBeats\\sample_data\\');

% Load EEG channel labels
labels = readtable('C:\Users\Tracy\Downloads\ds003838\sub-032\eeg\sub-032_task-memory_electrodes.tsv',"FileType","text",'Delimiter','\t');
for i = 1:size(labels,1)
    EEG.chanlocs(i).labels = labels.name{i};
end
locPath = fileparts(which('dipfitdefs.m'));
EEG = pop_chanedit(EEG,'lookup',fullfile(locPath,'standard_BEM','elec','standard_1005.elc'));

% Add ECG and PPG data (IMPORTANT: note that they were collected simultaneously
% with the EEG signals with the same hardware, avoiding potnetial time 
% synchronization issues). 
params.heart_channels = {CARDIO.chanlocs.labels};
EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
for iChan = 1:CARDIO.nbchan
    EEG.chanlocs(end+1).labels = params.heart_channels{iChan};
end
EEG = eeg_checkset(EEG);

% Downsample
EEG = pop_resample(EEG,250);

% Save
EEG = pop_saveset(EEG, 'filename','dataset1.set','filepath',fullfile(mainDir,'sample_data'));

% plot raw signals
pop_eegplot(EEG,1,1,1);
