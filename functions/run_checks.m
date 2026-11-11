% RUN_CHECKS - Check data and parameters, install missing plugins, set defaults.
%
% Usage:
%   [EEG, params, err] = run_checks(EEG, params, isCmd)
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

function [EEG, params, err] = run_checks(EEG, params, isCmd)

fprintf('Running basic checks... \n')

err = false;

% Only continuous data are supported (HEP epochs are created internally)
if length(size(EEG.data)) ~= 2
    bb_fail(isCmd, "Epoched EEG data detected. BrainBeats only supports continuous data at the moment.")
    err = true; return
end

% Sampling rate: at most 1024 Hz, 512 Hz recommended. Higher rates only cost
% memory and time (surrogate HEPs, ICA, nonlinear features); at 512 Hz the
% R-peaks are located within 2 ms. Above 1024 Hz: an error on the command
% line, a choice in the GUI; above 512 Hz: a warning, or a choice in the GUI.
if EEG.srate > 512
    [EEG, params, err] = check_srate(EEG, params, isCmd);
    if err, return; end
end

% HEP without EEG cleaning: the data must be high-pass filtered. The HEP is
% not baseline-corrected, so each channel's offset would be in the HEP, its
% scalp maps and the QRS map. Offset ratio: median |channel mean| / median
% channel SD within 1 s. On 7 recordings (2026-11-05): raw 105-1080;
% high-passed at 0.5 Hz below 0.6 (mostly < 0.003); at 0.1 Hz up to 9 (edge
% transients of large artifacts). Threshold 30, between the two.
if isfield(params,'analysis') && strcmp(params.analysis,'hep') && isfield(params,'clean_eeg') && ~params.clean_eeg
    r = offset_ratio(EEG, params);
    if r > 30
        bb_fail(isCmd, sprintf(['The EEG data do not look high-pass filtered: the channels'' offsets are %.3g times ' ...
            'their variability within 1 s (high-passed data: below 10; raw data: above 100). The HEP is not baseline-corrected, so these ' ...
            'offsets would be in the HEP, its scalp maps and the QRS map. High-pass filter the EEG first, e.g. ' ...
            'EEG = pop_eegfiltnew(EEG,''locutoff'',0.5); or set ''clean_eeg'' on (0.5-Hz high-pass for HEP).'], r));
        err = true; return
    end
end

% Check heart signal type
if ~contains(params.heart_signal, {'ecg' 'ppg' 'rr' 'off'})
    bb_fail(isCmd, 'Heart signal should be either ECG, PPG, RR, or off')
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
        bb_fail(isCmd, "At least one heart channel was not found in this dataset's channel list. Please make sure that you typed the correct label for your heart channels.")
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
        bb_fail(isCmd, "Electrode location coordinates must be loaded for visualizing outputs.")
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
    % GEDAI artifact removal (external plugin, noncommercial license)
    % asks first in the MATLAB desktop (license); installs directly without
    % one (e.g. -batch), where 'clean_method','gedai' is the request
    if isfield(params,'clean_method') && strcmpi(params.clean_method,'gedai') && ~exist('GEDAI','file')
        msg = ['GEDAI is not installed. It is an EEGLAB plugin under the PolyForm Noncommercial ' ...
            'license (no commercial use).'];
        fprintf('%s Installing it from the EEGLAB extension manager... \n', msg)
        if usejava('desktop') && ~strcmp(questdlg([msg ' Install it now?'], 'GEDAI', 'Install', 'Cancel', 'Install'), 'Install')
            error("GEDAI not installed: use 'clean_method','asr_ica', or install GEDAI from EEGLAB (File > Manage extensions).")
        end
        plugin_askinstall('GEDAI', 'GEDAI', 1);
        if ~exist('GEDAI','file')
            error("The GEDAI plugin could not be installed: install it from EEGLAB (File > Manage extensions), or use 'clean_method','asr_ica'.")
        end
    end
end

% Fuzzy entropy and fractal dimension (nonlinear HRV and EEG features) come
% from the ASCENT plugin. If it cannot be installed, BB_FUZZY_FRACTAL warns
% and sets them to NaN; the other features are computed.
if (isfield(params,'hrv_nonlinear') && params.hrv_nonlinear) || ...
        (isfield(params,'eeg_nonlinear') && params.eeg_nonlinear)
    if ~exist('compute_FuzzEn','file')
        try
            plugin_askinstall('ASCENT', 'compute_FuzzEn', 1);
        catch ME
            fprintf('The ASCENT plugin could not be installed (%s).\n', ME.message)
        end
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
        % ahead of the folders added on the client (the ASCENT plugin's
        % compute_psd took over on the workers: hence bb_compute_psd).
        % ASCENT's functions after it (a pool opened before ASCENT was
        % installed would not have them)
        pool = gcp('nocreate');
        if ~isempty(pool)
            bbpath = fileparts(which('eegplugin_BrainBeats.m'));
            folders = {fullfile(bbpath,'functions'), bbpath};
            if exist('compute_FuzzEn','file')
                folders{end+1} = fileparts(which('compute_FuzzEn'));
            end
            wait(parfevalOnAll(pool, @addpath, 0, folders{:}));
        end
    else
        warning("You do not have the parallel toolbox installed. Turning parallel computing OFF.")
        bb_warndlg("You do not have the parallel toolbox installed. Turning parallel computing OFF.")
        params.parpool = false;
    end
else
    fprintf('Parallel computing set to OFF. \n')
end

if ~isfield(params,'gong')
    params.gong = true;
end


function [EEG, params, err] = check_srate(EEG, params, isCmd)
% Ask to downsample above 512 Hz (GUI), or error above 1024 Hz and warn
% above 512 Hz (command line). Under 'matlab -batch' the GUI answer is the
% recommended one (512 Hz above 1024 Hz, else keep).
err = false;
fs = EEG.srate;
msg = sprintf(['The data are sampled at %g Hz. BrainBeats needs at most 1024 Hz, and 512 Hz is recommended: ' ...
    'higher rates only cost memory and time (at 512 Hz, R-peaks are located within 2 ms).'], fs);
batch = exist('batchStartupOptionUsed') > 0 && batchStartupOptionUsed; %#ok<EXIST>
if isCmd
    if fs > 1024
        bb_fail(isCmd, sprintf('%s Downsample first, e.g. EEG = pop_resample(EEG, 512);', msg));
    end
    warning('BrainBeats:srate', '%s To downsample: EEG = pop_resample(EEG, 512);', msg)
    return
end
rec = '512 Hz (recommended)';
if batch
    answer = rec; if fs <= 1024, answer = 'Keep'; end
    fprintf('%s Answer: %s. \n', msg, answer);
elseif fs > 1024
    answer = questdlg([msg ' Downsample to:'], 'Sampling rate', rec, '1024 Hz', 'Cancel', rec);
else
    answer = questdlg([msg ' Downsample to 512 Hz?'], 'Sampling rate', rec, sprintf('Keep %g Hz', fs), 'Cancel', rec);
end
switch strtok(answer)
    case {'512' '1024'}
        target = str2double(strtok(answer));
        fprintf('Downsampling from %g to %g Hz... \n', fs, target);
        EEG = pop_resample(EEG, target);
        params.resampled = [fs target];    % Hz: from, to (EEGLAB history)
    case 'Keep'
    otherwise
        disp('Aborted.'); err = true;
end


function r = offset_ratio(EEG, params)
% Median over EEG channels of |channel mean|, relative to the median over
% channels of their median SD within 1-s windows (heart and auxiliary
% channels excluded)
lab = {EEG.chanlocs.labels};
eeg = ~contains(lower(lab), {'ecg' 'ekg' 'ppg' 'aux' 'gsr' 'eda' 'eog' 'emg'});
if isfield(params,'heart_channels') && ~isempty(params.heart_channels)
    eeg = eeg & ~ismember(lab, cellstr(params.heart_channels));
end
x = double(EEG.data(eeg,:));
w = max(1, round(EEG.srate)); nW = floor(size(x,2)/w);
if isempty(x) || nW < 1, r = 0; return; end
sdw = median(std(reshape(x(:,1:nW*w), size(x,1), w, nW), 0, 2), 3);
r = median(abs(mean(x,2))) / max(median(sdw), eps);
