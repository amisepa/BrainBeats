% BrainBeats tutorial
% Launch each section one by one by clicking in the section and pressing: 
% CTRL/CMD + ENTER

clear; close all; clc
eeglab; close; % launch EEGLAB wihtout the GUI

% Find path to the plugin to load sample data files by command line 
mainDir = fileparts(which('eegplugin_BrainBeats.m')); 
cd(mainDir); % go to the plugin directory


%% METHOD 1: Heartbeat-evoked potentials (HEP) or oscillations (HEO) 

% Load sample raw data file into EEGLAB
EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true,'save',false,'vis',true); 

% Process file for HEP analysis using default parameters (WARNINGL EEG data 
% must be preprocessed because they will not be preprocessed by default!!!)
EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'}); 
% pop_eegplot(EEG,1,1,1);  % to visualize

% Same but adjusting some paramters this time
%   - 'clean_rr' to choose the method to interpolate RR artifacts (default = 'pchip')
%   - 'clean_eeg' to clean the EEG data (true) or not (false) 
%   - 'parpool' to turn ON (true) or OFF (false) the parrallel toolbox
%   - 'gpu' to turn ON (true) or OFF (false) GPU computing
%   - 'save' to save the resulting file (true) or not (false)
%   - 'vis' to visualize outputs with plots (true) or not (false)
EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'},'clean_rr','spline','clean_eeg',true, ...
    'parpool',false,'gpu',true,'save',false,'vis',true); 


%% METHOD 2: Extract EEG and HRV features

EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG1' 'ECG2'}, 'clean_rr','pchip','clean_eeg',true,'norm',true,...
    'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
    'eeg_features', {'time' 'frequency'}, ...
    'gpu',false,'parpool',true,'save',false,'vis',true);

%% METHOD 3: Remove heart components from EEG signals

EEG = pop_loadset('filename','sample_data2.set','filepath',fullfile(mainDir,'sample_data'));
EEG = brainbeats_process(EEG,'analysis','rm_heart','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'save',false,'vis',true);

%% To launch GUI via command line

EEG = brainbeats_process(EEG);

%% GROUP ANALYSIS: HEP

clear; close all; clc
eeglab; close;

% Download the open-source dataset here (4.3 GB): 
% https://nemar.org/dataexplorer/detail?dataset_id=ds001787
% Unzip the folder in the destination of your choice, and add the path to the
% directory below. Find details on the study here: 
% https://cerco.cnrs.fr/wp-content/uploads/2019/04/brandmeyer_t_18_2519.pdf
% data_dir = 'copy-paste you path here';
data_dir = 'C:\Users\Tracy\Documents\MATLAB\data_meditation_tracy';
cd(data_dir)

% Import the whole study using pop_importbids
[STUDY, ALLEEG] = pop_importbids( data_dir, 'eventtype','value', ...
    'bidsevent','on','bidschanloc','on', ...
    'outputdir',fullfile(data_dir,'derivatives') );
[STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);
EEG = ALLEEG; CURRENTSTUDY = 1; CURRENTSET = 1:length(EEG);

% Add session variable
STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name','STUDY.design 1','delfiles','off', ...
    'defaultdesign','off','variable1','group','values1',{'expert','novice'}, ...
    'vartype1','categorical','variable2','session','values2',{1,2,3},'vartype2', ...
    'continuous','subjselect',{ALLEEG.subject});

% Extract relevant info from STUDY
ids = {ALLEEG.subject};
% subjects = dir;
% subjects = {subjects.name}';
% subjects(~contains(subjects,'sub')) = [];
nSub = length(ids);  % n
filepaths = {ALLEEG.filepath};
filenames = {ALLEEG.filename};


design = STUDY.design;

% Process files for HEP analysis
for iSub = 1:nSub
    
    % Get subject data
    % EEG = pop_biosig(fullfile(filepath, filename));


    % Probes 
    % Q1: “Please rate the depth of your meditation,” for
    %   which participants evaluated the subjective depth of their
    %   “meditative state” for the moments immediately preceding
    %   the first probe, on a scale from 0 (not meditating at all) to
    %   3 (deep meditative state) by pressing the corresponding key
    %   on the keypad. 
    % Q2: “Please rate the depth of your mind wandering” automatically 
    %   followed. Participants evaluated the subjective depth of their 
    %   “mind wandering” for the period
    %   of time which immediately preceded the first probe, on a
    %   scale from 0 (not mind wandering at all) to 3 (immersed
    %   in their thoughts). 
    % Q3: "Please rate how tired you are,” where participants were asked to
    %   rate the subjective depth of their drowsiness at the time of the
    %   first question, from 0 (not drowsy at all) to 3 (very drowsy).

    % Markers
    % "2": "Response 1 (this may be a response to question 1, 2 or 3)",
    % "4": "Response 2 (this may be a response to question 1, 2 or 3)",
    % "8": "Response 3 (this may be a response to question 1, 2 or 3)",
    % "16": "Indicate involuntary response",
    % "128": "First question onset (most important marker)"

    % Load file
    filepath = fullfile(data_dir, subjects{iSub}, 'ses-01','eeg');
    % filename = sprintf('%s', subjects{iSub}, 'ses-01', );
    % sub-001_ses-01_task-meditation_eeg.bdf


    % Data were then segmented into 10 s-epochs from −10.05 to −0.05 s 
    % prior to onset of Q1 in the experience-sampling probe series.
    
    
    % You will be prompted to install the biosig plugin if not already installed
    % to load the .bdf files


end












