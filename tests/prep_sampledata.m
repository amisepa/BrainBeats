% PREP_SAMPLEDATA - Script that built the BrainBeats sample dataset
% (sample_data/dataset.set).
%
% Merges one subject's resting-state EEG with its ECG and PPG channels,
% adds 10-20 channel locations (standard_1005.elc), downsamples to 125 Hz,
% and adds artifacts so that the cleaning steps have something to remove:
% a bad channel (channel 10 time-reversed and x3), a large artifact on all
% EEG channels over the first 150 samples, and 3 s of white noise (SD 100)
% on temporal channels starting at 10 s to mimic muscle artifacts.
%
% Kept as a record of how the sample data were made: the hardcoded input
% and output paths are those of the machine that ran it (output saved as
% dataset-new.set; the repository copy is sample_data/dataset.set).
%
% Copyright (C) - Cedric Cannard, 2023

clear; close all; clc
datapath = 'C:\Users\Tracy\Downloads';
eeglab;close
locPath = fileparts(which('dipfitdefs.m'));
cd(datapath)

subject = 'sub-032';   % other available subjects: 32-36, 38-65

% Load EEG file
EEG = pop_loadset('filename',sprintf('%s_task-rest_eeg.set',subject),'filepath',datapath);

% Load 10-20 channel locations
EEG = pop_chanedit(EEG,'lookup',fullfile(locPath,'standard_BEM','elec','standard_1005.elc'));

% Load file containing ECG and PPG data
CARDIO = pop_loadset('filename',sprintf('%s_task-rest_ecg.set',subject),'filepath',datapath);

% if CARDIO data are a few seconds longer, select only common signal with EEG
extraData = CARDIO.xmax - EEG.xmax;
if extraData ~= 0
    warning("Removing %g seconds of extra data from cardio signal!",round(extraData,2))
    CARDIO = pop_select(CARDIO,'point',[1 EEG.pnts]);
end

% Check time is the same for both files
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
EEG = pop_resample(EEG,125);

% Artificially create a bad EEG channel
EEG.data(10,:) = EEG.data(10,end:-1:1).*3;

% Simulate a large electrode disconnection artifact at the beginning of
% the file (EEG channels only)
EEG.data(1:EEG.nbchan-length(heart_channels),1:150) = EEG.data(1:EEG.nbchan-length(heart_channels),150:-1:1).*3;

% Simulate high-frequency muscle artifacts 
channels = [9 10 20 21 42 55];      % temporal channels
startTime = 10*EEG.srate;           % Start at 10 s
duration = 3;                       % lasts 3 s
t = 0:1/EEG.srate:duration-1/EEG.srate;
artifact = 100 .* randn(length(channels), length(t));
EEG.data(channels, startTime:(startTime+length(t)-1)) = EEG.data(channels, startTime:(startTime+length(t)-1)) + artifact;

pop_saveset(EEG, 'filename','dataset-new.set','filepath','C:\Users\Tracy\Documents\MATLAB\BrainBeats\sample_data\');

% Remove DC offset and plot for visual inspection (not saved)
EEG.data = EEG.data - mean(EEG.data,2);
pop_eegplot(EEG,1,1,1);
