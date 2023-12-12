%% prep sample data for BrainBeats tutorial

datapath = 'C:\Users\Tracy\Downloads\ds003838';
eeglab;close

subject = 'sub-034';   %32-36, 38-44
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

% artifically create a bad EEG channel
EEG.data(10,:) = EEG.data(10,end:-1:1).*3;

% pop_eegplot(EEG,1,1,1);

pop_saveset(EEG, 'filename','dataset1.set','filepath','C:\Users\Tracy\Documents\\MATLAB\BrainBeats\sample_data\');

%% Simulate heart artifacts for method 3 by averaging heart signal into EEG signals

EEG = pop_loadset('dataset1.set','C:\Users\Tracy\Documents\\MATLAB\BrainBeats\sample_data\');
EEG = pop_eegfiltnew(EEG,'locutoff',1);
EEG = pop_eegfiltnew(EEG,'hicutoff',30);
CARDIO = pop_select(EEG,'channel',{'ECG'}); 
EEG = pop_select(EEG,'nochannel',{'PPG' 'ECG'}); 
EEG = ref_infinity(EEG);

heart_channels = {CARDIO.chanlocs.labels};
for iChan = 1:CARDIO.nbchan
    CARDIO.data(iChan,:) = rescale(CARDIO.data(iChan,:), -500, 500);
end
EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
for iChan = 1:CARDIO.nbchan
    EEG.chanlocs(end+1).labels = heart_channels{iChan};
end
EEG = eeg_checkset(EEG);

EEG.data = EEG.data - repmat(mean(EEG.data),size(EEG.data,1),1);
pop_eegplot(EEG,1,1,1);

EEG = pop_saveset(EEG, 'filename','dataset2.set','filepath','C:\Users\Tracy\Documents\MATLAB\BrainBeats\sample_data\');
