%% Get parameters specified by user via command line
% 
% Cedric Cannard, April 2023

function params = getparams_command(varargin)

inputs = varargin(1:2:end);

% Extract user parameters
idx = find(contains(inputs,'heart_signal'));
if ~isempty(idx)
    params.heart_signal = lower(varargin{idx*2});
    if ~contains(params.heart_signal, {'ecg' 'ppg'})
        error("Heart signal not recognized. Should be 'ECG' or 'PPG'.")
    end
else
    error("Heart signal type not defined. Please define 'heart_signal' as 'ECG' or 'PPG'. See help for an example")
end

% Detect if there are 1 or multiple heart channel(s)
idx = find(contains(inputs,'heart_channels'));
if ~isempty(idx)
    params.heart_channels = varargin{idx*2};
else
    error("Heart channels not defined. Please define 'heart_channels'. See help for an example")
end

% Analysis to do
idx = find(contains(inputs,'analysis'));
if ~isempty(idx)
    params.analysis = varargin{idx*2};
end

% RR artifact correction method
idx = find(contains(inputs,'rr_correct'));
if ~isempty(idx)
    params.rr_correct = varargin{idx*2};
elseif isempty(idx) && contains(params.analysis, {'features' 'hep'})
    fprintf('RR artifact correction method not defined. Selecting Linear interpolation by default. \n')
    params.rr_correct = 'linear';
end

% Clean EEG
idx = find(contains(inputs,'clean_eeg'));
if ~isempty(idx)
    params.clean_eeg = varargin{idx*2};
else
    params.clean_eeg = false;
end

% EEG features
if strcmp(params.analysis,'features')
    idx = find(contains(inputs,'eeg_features'));
    if ~isempty(idx)
        params.eeg = true;
        eeg_features = varargin{idx*2};
        if sum(contains(eeg_features,{'frequency'})) > 0
            params.eeg_frequency = true;
        else
            params.eeg_frequency = false;
        end
        if sum(contains(eeg_features,'nonlinear')) > 0
            params.eeg_nonlinear = true;
        else
            params.eeg_nonlinear = false;
        end            
    else
        params.eeg = false;
    end
end

% HRV features
if strcmp(params.analysis,'features')
    idx = find(contains(inputs,'hrv_features'));
    if ~isempty(idx)
        params.hrv = true;
        hrv_features = varargin{idx*2};
        if sum(contains(hrv_features,{'time'})) > 0
            params.hrv_time = true;
        else
            params.hrv_time = false;
        end
        if sum(contains(hrv_features,{'frequency'})) > 0
            params.hrv_frequency = true;
        else
            params.hrv_frequency = false;
        end
        if sum(contains(hrv_features,'nonlinear')) > 0
            params.hrv_nonlinear = true;
        else
            params.hrv_nonlinear = false;
        end            
    else
        params.hrv = true;
    end
else
    params.hrv = false;
end

% Data visualization
idx = find(contains(inputs,'vis'));
if ~isempty(idx)
    params.vis = logical(varargin{idx*2});
else
    disp('Visualization not defined. Visualization is turned ON by default')
    params.vis = true;
end

% GPU computing
idx = find(contains(inputs,'gpu'));
if ~isempty(idx)
    params.gpu = logical(varargin{idx*2});
else
    params.gpu = false;
end
