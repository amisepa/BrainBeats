%% prep sample data for BrainBeats tutorial

datapath = 'C:\Users\Tracy\Downloads\ds003838';
eeglab;close

subject = 'sub-032';   %32-36, 38-44
datapath = fullfile(datapath,subject);
cd(datapath)

% Load EEG file
EEG = pop_loadset('filename',sprintf('%s_task-rest_eeg.set',subject),'filepath',fullfile(datapath,'eeg'));

% Load 10-20 channel locations
locPath = fileparts(which('dipfitdefs.m'));
EEG = pop_chanedit(EEG,'lookup',fullfile(locPath,'standard_BEM','elec','standard_1005.elc'));

% Load ECG and PPG file
CARDIO = pop_loadset('filename',sprintf('%s_task-rest_ecg.set',subject),'filepath',fullfile(datapath,'ecg'));

% Pull cardio channel names
heart_channels = {CARDIO.chanlocs.labels};

% Merge
if EEG.pnts ~= CARDIO.pnts
    EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data(:,1:end-1); % for PPG
else
    EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
end
EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
for iChan = 1:CARDIO.nbchan
    EEG.chanlocs(end+1).labels = heart_channels{iChan};
end
EEG = eeg_checkset(EEG);

% downsample so that the repo is not too heavy and computations are fast
EEG = pop_resample(EEG,250);

% Artifically create a bad EEG channel
EEG.data(10,:) = EEG.data(10,end:-1:1).*3;

% Simulate a large electrode disconnection artifact in the beginning of
% file
EEG.data(1:63,1:150) = EEG.data(1:63,150:-1:1).*3;

% Simulate high-frequency muscle artifacts 
channels = [9 10 20 21 42 55];  % TP channels
startTime = 10*EEG.srate;        % Start at 3 s
duration = 3;                   % duration of artifact (in s)
t = 0:1/EEG.srate:duration-1/EEG.srate;
artifact = 100 .* randn(length(channels), length(t));
EEG.data(channels, startTime:(startTime+length(t)-1)) = EEG.data(channels, startTime:(startTime+length(t)-1)) + artifact;

pop_saveset(EEG, 'filename','dataset.set','filepath','C:\Users\Tracy\Documents\\MATLAB\BrainBeats\sample_data\');

% EEG = rm_DC(EEG);
% pop_eegplot(EEG,1,1,1);

%% Simulate heart artifacts for method 3 by averaging heart signal into EEG signals

% EEG = pop_loadset('dataset.set','C:\Users\Tracy\Documents\\MATLAB\BrainBeats\sample_data\');
% EEG = pop_eegfiltnew(EEG,'locutoff',1);
% EEG = pop_eegfiltnew(EEG,'hicutoff',30);
% CARDIO = pop_select(EEG,'channel',{'ECG'}); 
% EEG = pop_select(EEG,'nochannel',{'PPG' 'ECG'}); 
% EEG = ref_infinity(EEG);
% 
% heart_channels = {CARDIO.chanlocs.labels};
% for iChan = 1:CARDIO.nbchan
%     CARDIO.data(iChan,:) = rescale(CARDIO.data(iChan,:), -500, 500);
% end
% EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
% EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
% for iChan = 1:CARDIO.nbchan
%     EEG.chanlocs(end+1).labels = heart_channels{iChan};
% end
% EEG = eeg_checkset(EEG);
% 
% EEG.data = EEG.data - repmat(mean(EEG.data),size(EEG.data,1),1);
% pop_eegplot(EEG,1,1,1);
% 
% EEG = pop_saveset(EEG, 'filename','dataset2.set','filepath','C:\Users\Tracy\Documents\MATLAB\BrainBeats\sample_data\');
