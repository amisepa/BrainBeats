%% Welcome to the BrainBeats tutorial for command line use
% 
% REQUIREMENTS:
%   1) Install EEGLAB
%   Download here: https://github.com/sccn/eeglab
%   Unzip (or clone) the file on your computer and add the path to MATLAB 
%   Home panel > Set path > Add folder > select the eeglab folder > Save >
%   Close
% 
%   2) Install BrainBeats
%   Download here: https://github.com/amisepa/BrainBeats
%   Type 'eeglab' in MATLAB's command window to open EEGLAB. Go to File >
%   Manage extensions > type 'brainbeats' in the search bar > select in the
%   list area, and click Install. Or, if you use Git, simply clone the repo 
%   in eeglab > plugins on your computer. 
% 
% We hope you enjoy the tutorial! 
% 
% You can launch each section one by one by clicking in the section and pressing:
% CTRL/CMD + ENTER

%% Open EEGLAB and get the path to the plugin automatically

clear; close all; clc

% Launch EEGLAB without the GUI
eeglab; close;

% Get the path to the EEGLAB plugin which contains a folder with the sample data
main_path = fileparts(which('eegplugin_BrainBeats.m'));

% Go to the plugin directory
cd(main_path)

%% METHOD 1: Heartbeat-evoked potentials (HEP) and oscillations (HEO)

% Load the dataset1 into EEGLAB. This file contains raw 64-channel EEG,
% ECG, and PPG data during 3.8 minutes of resting state eyes opened. This
% file corresponds to sub-032_task-rest_eeg.set and
% sub-032_task-rest_ecg.set merged, available here:
% https://nemar.org/dataexplorer/detail?dataset_id=ds003838
% note: sub-032_task-rest_ecg.set contains both ECG and PPG channels
% note: we artificially modified channel 37 to be detected as a bad channel
% for demonstration purposes since there were no bad channels in this
% dataset. 
EEG = pop_loadset('filename','dataset1.set','filepath',fullfile(main_path,'sample_data'));

% Process file for HEP analysis using default parameters except for:
%   - selecting the type of analysis: 'hep'
%   - selecting the type of heart signal: 'ECG'
%   - selecting the name of the ECG electrodes: 'ECG'
%   - 'clean_eeg' set to to 'true' preprocess the EEG data with default
%       parameters.
% Note: the toolbox automatically detects the undesired PPG channel, 
% which is expected since the toolbox is not designed to run both ECG and PPG at the time.
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true);

%% Same as above but using the PPG signal and adjusting some parameters 
%  Note that we are changing these parameters for illustraiton only, but
%  default parameters are recommended. These parameters should only be
%  changed if you know why. 

% We need to load the file again since it was modified by our previous call above. 
EEG = pop_loadset('filename','dataset1.set','filepath',fullfile(main_path,'sample_data'));

% Here we change the following parameters:
%   - 'heart_signal' set to 'PPG' (signal type)
%   - 'heart_channels' set to 'PPG' (electrode name)
%   - 'clean_rr' set to 'pchip' to interpolate the RR artifacts
%   - 'reref' set to 'infinity' to rereference EEG data to infinity
%   - 'highpass' set to 2 to remove EEG frequencies <2 hz
%   - 'lowpass' set to 30 to remove EEG frequencies >80 hz
%   - 'filttype' set to noncausal to use noncausal zero-phase filter
%       instead of the default causal minimum-phase filter
%   - 'detectMethod' set to 'mean' to remove bad EEG epochs more
%       conservatively (default = 'grubbs')
%   - 'save' set to 'true' to save the final 'filename_HEP.set' file
%   - 'vis_cleaning' set to 'true' as we already generated all the preprocessing
%       plots above.
%   - 'vis_outputs' set to 'true' to visualize the final output again. 
% Note: the toolbox automatically detects the undesired ECG channel, 
% which is expected since the toolbox is not designed to run both ECG and PPG at the time.
EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'clean_rr','pchip','clean_eeg',true, ...
    'reref','infinity','highpass',2,'lowpass',30,'filttype','noncausal', ...
    'detectMethod','mean','save',true,'vis_cleaning',true,'vis_outputs',true);

%% METHOD 2: Extract EEG and HRV features

% Load the same raw dataset again
EEG = pop_loadset('filename','dataset1.set','filepath',fullfile(main_path,'sample_data'));

% This time, note these different inputs:
%   - 'analysis' set to 'features' to extract EEG and HRV features
%   - 'norm' to normalize the frequency domain outputs (type
%       'help get_hrv_features' and 'help get_eeg_features' for more detail)
%   - 'hrv_features' and 'eeg_features' to {'time' 'frequency' 'nonlinear'} 
%       to extract features in all three domains. (type 
%       'help get_hrv_features' and 'help get_eeg_features' for more detail)
%   - 'parpool' to 'true' to use parallel computing.
%   - 'vis_cleaning' to 'true' to see preprocessing steps (slightly different
%       than for HEP). 
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
    'eeg_features', {'time' 'frequency' 'nonlinear'}, ...
    'eeg_norm',1,'hrv_norm',false,'asy_norm',false, ...
    'gpu',false,'parpool',true,'save',false,...
    'vis_cleaning',true,'vis_outputs',true);

% You can find all features in EEG.brainbeats.features
% You can plot them again using:
% params.chanlocs = EEG.chanlocs;
% plot_features(EEG.brainbeats.features,params)

%% METHOD 3: Remove heart components from EEG signals

EEG = pop_loadset('filename','dataset1.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_eegfiltnew(EEG,'locutoff',1);
EEG = pop_eegfiltnew(EEG,'hicutoff',30);
EEG = brainbeats_process(EEG,'analysis','rm_heart','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'keep_heart',true, ...
    'save',true,'vis_outputs',true);

%% To launch the GUI only

EEG = pop_loadset('filename','dataset1.set','filepath',fullfile(main_path,'sample_data'));
EEG = brainbeats_process(EEG);

% Type 'eegh' at the end of the operations to output the command line with
% the parameters that were selected manually in the GUI
eegh

%% To process only cardiovascular signals and extract HRV features (no EEG)
% To turn OFF all EEG operations, the input 'eeg' is set to 'false'.

% ECG
EEG = pop_loadset('filename','dataset1.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'PPG'}); 
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'eeg',false,...
    'hrv_features',{'time' 'frequency' 'nonlinear'},'norm',false, ...
    'vis_cleaning',true,'vis_outputs',true);

% PPG
EEG = pop_loadset('filename','dataset1.set','filepath',fullfile(main_path,'sample_data'));
EEG = pop_select(EEG,'nochannel',{'ECG'}); 
EEG = brainbeats_process(EEG,'analysis','features','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'eeg',false,...
    'hrv_features',{'time' 'frequency' 'nonlinear'},...
    'vis_cleaning',true,'vis_outputs',true);

% Cardiovascular preprocessing outputs can be found in:
EEG.brainbeats.preprocessings

% Cardiovascular preprocessing outputs can be found in:
EEG.brainbeats.features.HRV

%% Group analysis: HEP



%% Group analysis: Features



%% Save UI figures in high definition

cd("figures")
figname = inputdlg('Figure name for saving:');

% UI figures
exportgraphics(gca,sprintf('%s.png',figname{:}),'resolution',300)
exportgraphics(gca,sprintf('%s.tiff',figname{:}),'resolution',300)

%% Save UI figures in high definition

figname = inputdlg('Figure name for saving:');

saveas(gcf,sprintf('%s.fig',figname{:}))
print(gcf, sprintf('%s.png',figname{:}),'-dpng','-r300');   % 300 dpi .png
print(gcf, sprintf('%s.tiff',figname{:}),'-dtiff','-r300');  % 300 dpi .tiff
