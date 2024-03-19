%% Run some basic checks to ensure basic parameters and requirements are
% met, install required plugins that are not already installed, and apply
% some default parameters when missing.
%
% Copyright (C) BrainBeats - Cedric Cannard 2023

function [EEG, params, err] = run_checks(EEG, params)

fprintf('Running basic checks... \n')

err = false;

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

% % Any heart operations?
% if ~isfield(params,'heart') || ~isempty(params.heart_signal) || ~isempty(params.heart_channels) || ...
%         isfield(params,'hrv_features') || strcmpi(params.analysis,'hep')
%     params.heart = true;
% else
%     params.heart = false;
% end

% Heart checks
% Make sure Heart channel is a cell
if ~iscell(params.heart_channels)
    % warning("Heart channel label should be a cell (e.g. {'ECG'} or {'AUX1' 'AUX2'}). Converting it to cell now.")
    params.heart_channels = {params.heart_channels};
end

% Check if heart channels are in file (for command line mode)
nchan = length(params.heart_channels);
idx = nan(nchan,1);
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

% Includes HRV or not (for plotting only)
if ~isfield(params,'hrv_features')
    % if strcmp(params.analysis,{'features'}) && params.clean_heart
        params.hrv_features = false;
    % end
end

% EEG checks
% if isfield(params, 'eeg')
% Includes EEG or not (for plotting only)
if ~isfield(params,'eeg_features')
    % if any(strcmp(params.analysis,{'hep' 'rm_heart'})) || (isfield(params, 'eeg_frequency') ...
    %         && params.eeg_frequency) || (isfield(params, 'eeg_nonlinear') && params.eeg_nonlinear)
    %     params.eeg_features = true;
    % else
        params.eeg_features = false;
        % params.clean_eeg = false;
    % end
end

% Check for channel locations
if params.eeg_features~=0
    if ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(2).X)
        errordlg("Electrode location coordinates must be loaded for visualizing outputs.")
        err = true; return
    end
end

% Set default ICA method if not already set
% 1 = picard (fast); 2 = infomax (default); 3 = modified infomax for replicability (long)
if ~isfield(params,'icamethod')
    params.icamethod = 2;
end

% EEG re-referencing
if ~isfield(params,'ref')
    params.ref = 'average';
end

% Install necessary plugins for preprocessing
if params.clean_eeg

    if ~exist('clean_asr','file')
        plugin_askinstall('clean_asr','clean_asr', 1);
    end
    if ~exist('picard','file') && params.icamethod == 1
        plugin_askinstall('picard', 'picard', 1);
    end
    if ~exist('iclabel','file')
        plugin_askinstall('iclabel', 'iclabel', 1);
    end
    if strcmp(params.ref, 'infinity')
        if ~exist('ref_infinity','file')
            plugin_askinstall('REST_cmd', 'REST_cmd', 1);
        end
    end
end
if strcmp(params.analysis,'rm_heart')
    if ~exist('picard','file') && params.icamethod == 1
        plugin_askinstall('picard', 'picard', 1);
    end
    if ~exist('iclabel','file')
        plugin_askinstall('iclabel', 'iclabel', 1);
    end
end
% else
% params.eeg = [];
% end

% Ensure data have double precision
EEG.data = double(EEG.data);

% Store sampling frequency
params.fs = EEG.srate;


% Initiate or block parallel computing
if params.parpool

    % check if user has parallel toolbox
    try
        ver('parallel')
    catch
        warning("You do not have the parallel toolbox. Turning parallel computing OFF.")
        warningdlg("You do not have the parallel toolbox. Turning parallel computing OFF.")
        params.parpool = false;
    end

    ps = parallel.Settings;
    fprintf('Parallel computing set to ON. \n')
    ps.Pool.AutoCreate = true;
    p = gcp('nocreate');
    % delete(gcp('nocreate')) % shut down opened parpool
    if isempty(p) % if not already on, launch it
        disp('Initiating parrallel computing (all cores and threads -1)...')
        c = parcluster; % cluster profile
        % N = feature('numcores');          % only physical cores
        N = getenv('NUMBER_OF_PROCESSORS'); % all processor (cores + threads)
        if ischar(N), N = str2double(N); end
        c.NumWorkers = N-2;  % update cluster profile to include all workers
        c.parpool();
    end
else
    fprintf('Parallel computing set to OFF. \n')

    % Shut down parallel pool if already opened
    p = gcp('nocreate');
    if ~isempty(p)
        delete(gcp('nocreate'));
    end

    % Prevent parfor loops from launching parpool mode
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false;  
end

if ~isfield(params,'gong')
    params.gong = true;
end
