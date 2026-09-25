% RUN_CHECKS - Check data and parameters, install missing plugins, set defaults.
%
% Usage:
%   [EEG, params, err] = run_checks(EEG, params)
%
% Checks (error dialog and err = true on failure): data are continuous;
% params.heart_signal is 'ecg', 'ppg', 'rr' or 'off'; all
% params.heart_channels exist in EEG.chanlocs (converted to a cell if
% needed); channel coordinates are present when EEG features are plotted.
% Offers to install clean_rawdata, ICLabel, Picard (icamethod 1) and
% REST_cmd (ref 'infinity') when needed. If params.parpool, starts a pool
% with the 'Processes' profile (cores - 1 workers) and puts BrainBeats first
% on the workers' path; turns params.parpool off if that fails.
%
% Inputs:
%   EEG    - EEGLAB EEG structure (EEG + heart channels)
%   params - BrainBeats parameters (reads heart_signal, heart_channels,
%            analysis, clean_eeg, vis_outputs, parpool)
% Outputs:
%   EEG    - with data converted to double
%   params - with fs = EEG.srate and defaults when missing: hrv_features =
%            false, eeg_features = false, icamethod = 2 (Infomax),
%            ref = 'average', gong = true
%   err    - true if a check failed (outputs then incomplete)
%
% Copyright (C) - Cedric Cannard, 2023

function [EEG, params, err] = run_checks(EEG, params)

fprintf('Running basic checks... \n')

err = false;

% Only continuous data are supported (HEP epochs are created internally)
if length(size(EEG.data)) ~= 2
    errordlg("Epoched EEG data detected. BrainBeats only supports continuous data at the moment.")
    err = true; return
end

% Check heart signal type
if ~contains(params.heart_signal, {'ecg' 'ppg' 'rr' 'off'})
    errordlg('Heart signal should be either ECG, PPG, RR, or off')
    err = true; return
end

% Heart checks
% Make sure Heart channel is a cell
if ~strcmpi(params.heart_signal,'off') && ~strcmpi(params.heart_signal,'rr')
    if ~iscell(params.heart_channels)
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
    end
end

% Includes HRV or not (for plotting only)
if ~isfield(params,'hrv_features')
    params.hrv_features = false;
end

% Includes EEG or not (for plotting only)
if ~isfield(params,'eeg_features')
    params.eeg_features = false;
end

% Check for channel locations
if params.eeg_features~=0 && params.vis_outputs
    if ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(1).X)
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
        plugin_askinstall('clean_rawdata','clean_asr', 1);
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

% Ensure data have double precision
EEG.data = double(EEG.data);

% Store sampling frequency
params.fs = EEG.srate;


% Initiate parallel computing if requested (parfor loops run serially
% otherwise: see get_eeg_features)
addons = ver;
parpool_installed = any(contains({addons.Name}, 'Parallel'));
if params.parpool
    if parpool_installed
        fprintf('Parallel computing set to ON. \n')
        if isempty(gcp('nocreate')) % if not already on, launch it
            N = feature('numcores');    % physical cores (works on all platforms)
            fprintf('Initiating parallel computing (%g workers)...\n', max(1,N-1))
            try
                % 'Processes' profile explicitly: the default profile may be
                % 'Threads', which parcluster does not accept
                c = parcluster('Processes');
                c.NumWorkers = max(1,N-1);
                parpool(c, c.NumWorkers);
            catch ME
                warning('Could not start a parallel pool (%s). Turning parallel computing OFF.', ME.message)
                params.parpool = false;
            end
        end

        % Put BrainBeats first on the workers' path too: they start from the
        % saved path, which can hold other functions with the same names
        % (e.g. another compute_psd) ahead of the folders added on the client
        pool = gcp('nocreate');
        if ~isempty(pool)
            bbpath = fileparts(which('eegplugin_BrainBeats.m'));
            wait(parfevalOnAll(pool, @addpath, 0, fullfile(bbpath,'functions'), bbpath));
        end
    else
        warning("You do not have the parallel toolbox installed. Turning parallel computing OFF.")
        warndlg("You do not have the parallel toolbox installed. Turning parallel computing OFF.")
        params.parpool = false;
    end
else
    fprintf('Parallel computing set to OFF. \n')
end

if ~isfield(params,'gong')
    params.gong = true;
end
