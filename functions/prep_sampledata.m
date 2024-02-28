%% prep sample data for BrainBeats tutorial

clear; close all; clc
datapath = 'C:\Users\Tracy\Downloads';
eeglab;close
locPath = fileparts(which('dipfitdefs.m'));
cd(datapath)

subject = 'sub-032';   %32-36, 38-65

% Load EEG file
EEG = pop_loadset('filename',sprintf('%s_task-rest_eeg.set',subject),'filepath',datapath);

% Load 10-20 channel locations
EEG = pop_chanedit(EEG,'lookup',fullfile(locPath,'standard_BEM','elec','standard_1005.elc'));

% Load file containing ECG and PPG data
CARDIO = pop_loadset('filename',sprintf('%s_task-rest_ecg.set',subject),'filepath',datapath);

% CARDIO data are a few seconds longer, select only common signal with EEG
extraData = CARDIO.xmax - EEG.xmax;
warning("Removing %g seconds of extra data from cardio signal!",round(extraData,2))
CARDIO = pop_select(CARDIO,'point',[1 EEG.pnts]);

% Check time is the same for both files
% CARDIO.times = double(CARDIO.times);
tmp = diff([single(EEG.times); single(CARDIO.times)]);
if any(tmp~=0)
    nSamples = tmp~=0;
    error("%g%% of the samples have a different time stamp between CARDIO and EEG data! Meaning time synchronization will not be ensured", round(nSamples/EEG.pnts*100,2))
end

% Pull cardio channel names
heart_channels = {CARDIO.chanlocs.labels};

% Merge
EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
for iChan = 1:CARDIO.nbchan
    EEG.chanlocs(end+1).labels = heart_channels{iChan};
end
EEG = eeg_checkset(EEG);

% downsample so that the repo is not too heavy and computations are fast
EEG = pop_resample(EEG,250);

% Artifically create a bad EEG channel
EEG.data(10,:) = EEG.data(10,end:-1:1).*5;

% Simulate a large electrode disconnection artifact in the beginning of
% file
EEG.data(1:EEG.nbchan-length(heart_channels),1:EEG.srate) = EEG.data(1:EEG.nbchan-length(heart_channels),EEG.srate:-1:1).*3;

% Simulate high-frequency muscle artifacts 
channels = [9 10 20 21 42 55];  % TP channels
startTime = 10*EEG.srate;        % Start at 10 s
duration = 3;                   % duration of artifact (in s)
t = 0:1/EEG.srate:duration-1/EEG.srate;
artifact = 100 .* randn(length(channels), length(t));
EEG.data(channels, startTime:(startTime+length(t)-1)) = EEG.data(channels, startTime:(startTime+length(t)-1)) + artifact;

pop_saveset(EEG, 'filename','dataset-new.set','filepath','C:\Users\Tracy\Documents\\MATLAB\BrainBeats\sample_data\');

% EEG = rm_DC(EEG);
pop_eegplot(EEG,1,1,1);

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
