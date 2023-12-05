%% Run some basic checks to ensure basic parameters and requirements are
% met, install required plugins that are not already installed, and apply
% some default parameters when missing.
%
% Copyright (C) BrainBeats - Cedric Cannard 2023

function [EEG, params, err] = run_checks(EEG, params)

fprintf('Running basic checks... \n')

err = false;

if isempty(EEG.ref)
    warning('EEG data not referenced! Referencing is highly recommended');
end

% Check if data format is compatible with chosen analysis and select analysis
% if isfield(params,'analysis')
% switch params.analysis
% case 'features'
if length(size(EEG.data)) ~= 2
    errordlg("Epoched EEG data detected. BrainBeats only supports continuous data at the moment.")
    err = true; return
end
% case 'epoched'
%     if length(size(EEG.data)) ~= 3
%         error("You selected HEP analysis but your data are not epoched.")
%     end
% end
% else
%     % Select analysis based on data format if not defined
%     if length(size(EEG.data)) == 2
%         params.analysis = 'continuous';
%         disp("Analysis not defined. Continuous data detected: selecting 'feature-based mode' by default")
%     elseif length(size(EEG.data)) == 3
%         params.analysis = 'epoched';
%         disp("Analysis not defined. Epoched data detected: selecting 'heart-beat evoked potential (HEP) mode' by default")
%     else
%         error("You did not define the analysis to run, and your data format was not recognized. " + ...
%             "Should be 'continuous' or 'epoched', and something may be wrong with your data format ")
%     end
% end

% Make sure Heart channel is a cell
if ~iscell(params.heart_channels)
    warning("Heart channel label should be a cell (e.g. {'ECG'} or {'AUX1' 'AUX2'}). Converting it to cell now.")
    params.heart_channels = {params.heart_channels};
end

% Check if heart channels are in file (for command line mode)
nchan = length(params.heart_channels);
for i = 1:nchan
    idx(i) = any(strcmp(params.heart_channels{i},{EEG.chanlocs.labels}));
    if idx(i) == 0
        warning("Heart channel %s not found in this dataset's channel list.",params.heart_channels{i})
    end
end
if length(idx) ~= sum(idx) 
    errordlg("At least one heart channel was not found in this dataset's channel list. Please make sure that you typed the correct label for your heart channels.")
    err = true; return
    % else
    %     fprintf("%g/%g heart channels confirmed in this dataset's channel list. \n", sum(idx), length(params.heart_channels))
% elseif length(idx) == 1 && sum(idx) == 0
%     errordlg("The heart channel label you typed was not found in this dataset's channel list. Please make sure that you typed the correct label for your heart channel.");
%     err = true; return
end

% Check heart signal type
if ~contains(params.heart_signal, {'ecg' 'ppg'})
    errordlg('Heart signal should be either ECG or PPG')
    return
end

% Includes EEG or not (for plotting only)
if ~isfield(params,'eeg')
    if any(strcmp(params.analysis,{'hep' 'rm_heart'})) || params.eeg_frequency
        params.eeg = true;
    else
        params.eeg = false;
        params.clean_eeg = false;
    end
end

% Check for channel locations
if params.eeg
    if ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(2).X)
        errordlg("Electrode location coordinates must be loaded for visualizing outputs.")
        err = true; return
    end
end

% Install necessary plugins for preprocessing
if params.clean_eeg
    % fprintf('Checking if required EEGLAB plugins are installed... \n')
    if ~exist('clean_asr','file')
        plugin_askinstall('clean_asr','clean_asr', 0);
    end
    if ~exist('picard','file')
        plugin_askinstall('picard', 'picard', 0);
    end
    if ~exist('iclabel','file')
        plugin_askinstall('iclabel', 'iclabel', 0);
    end
    if ~exist('ref_infinity','file')
        plugin_askinstall('REST_cmd', 'REST_cmd', 0);
    end
end
if strcmp(params.analysis,'rm_heart')
    if ~exist('iclabel','file')
        plugin_askinstall('iclabel', 'iclabel', 0);
    end
end


% Ensure data double precision
EEG.data = double(EEG.data);  % ensure double precision

% Store sampling frequency
params.fs = EEG.srate;

