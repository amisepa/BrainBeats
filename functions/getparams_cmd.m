% GETPARAMS_CMD - Parse the 'key', value inputs of BRAINBEATS_PROCESS into
% the params structure, check them, and set the defaults that depend on
% the analysis. The other defaults are set in RUN_CHECKS and in the
% functions that use them.
%
% Usage:
%   params = getparams_cmd('key', val, ...)
%
% See BRAINBEATS_PROCESS for the list of inputs.
%
% Copyright (C) - Cedric Cannard, 2023

function params = getparams_cmd(varargin)

%% Inputs

% 'key', value pairs: values are never read as keys ('ppg_transit','eeg' is
% not the 'eeg' input)
if mod(numel(varargin),2) ~= 0
    error('BrainBeats:input', "BrainBeats inputs must be 'key', value pairs.")
end
keys = varargin;
keys(2:2:end) = {''};
given = cellfun(@char, varargin(1:2:end), 'UniformOutput', false);
unknown = given(~ismember(lower(given), lower(known_keys())));
if ~isempty(unknown)
    error('BrainBeats:input', "Unknown input(s): %s. See 'help brainbeats_process'.", strjoin(unknown, ', '))
end
legacy = given(ismember(lower(given), legacy_keys()));
if ~isempty(legacy)
    warning('Input(s) ignored, no longer used: %s.', strjoin(legacy, ', '))
end

%% General parameters

% Extract user parameters
idx = find(strcmpi(keys,'heart_signal'));
if ~isempty(idx)
    params.heart = true;
    params.heart_signal = lower(varargin{idx+1});
    if ~contains(params.heart_signal, {'ecg' 'ppg' 'rr' 'off'})
        error("Heart signal not recognized. Should be 'ecg', 'ppg', 'rr', or 'off'.")
    end
    if contains(params.heart_signal, 'off')
        params.heart = false;
        params.heart_channels = {};
    end
else
    error("Heart signal type not defined. Please define 'heart_signal' as 'ecg', 'ppg', 'rr', or 'off'. Type help in the command window for an example")
end

% Detect if there are 1 or multiple heart channel(s)
if params.heart
    idx = find(strcmpi(keys,'heart_channels'));
    if ~isempty(idx)
        params.heart_channels = varargin{idx+1};
        fprintf('Number of heart channels selected: %g \n', length(params.heart_channels));
    elseif ~strcmpi(params.heart_signal, 'rr')
        error("Heart channels not defined. Please define 'heart_channels'. See help for an example")
    else
        params.heart_channels = {};
    end

    % Pre-detected beat latencies (for 'rr' mode)
    idx = find(strcmpi(keys,'beat_latencies'));
    if ~isempty(idx)
        params.beat_latencies = varargin{idx+1};
        fprintf('Pre-detected beat latencies provided: %g beats \n', length(params.beat_latencies));
    elseif strcmpi(params.heart_signal, 'rr')
        error("'beat_latencies' required when heart_signal is 'rr'. Provide a vector of beat times in seconds.")
    end
end

% Analysis to do
idx = find(strcmpi(keys,'analysis'));
if ~isempty(idx)
    params.analysis = lower(varargin{idx+1});
    fprintf('Selected mode to run: %s \n', params.analysis);
else
    error("Analysis type not defined. Must be set to either 'hep', 'features', or 'coherence'. ")
end
if strcmp(params.analysis, 'rm_heart')
    error('BrainBeats:input', ['The ''rm_heart'' analysis was removed in v1.6: ICA with the ECG channel ' ...
        'removed most of a synthetic heartbeat-locked brain response on 5 of 6 test subjects. Use ' ...
        '''analysis'',''hep'' with ''heart_removal'' (''ica'' or ''ecg_regression'').'])
end
if ~any(strcmp(params.analysis, {'hep' 'features' 'coherence'}))
    error("Analysis '%s' not recognized. Must be 'hep', 'features', or 'coherence'.", params.analysis)
end
if ~params.heart && any(strcmp(params.analysis, {'hep' 'coherence'}))
    error("The '%s' analysis needs a heart signal: 'heart_signal' cannot be 'off'.", params.analysis)
end
if strcmp(params.heart_signal,'rr') && strcmp(params.analysis, 'coherence')
    error("The '%s' analysis needs the ECG/PPG signal itself: 'heart_signal' cannot be 'rr'.", params.analysis)
end

%% Cardiovascular signal preprocessing

if params.heart
    idx = find(strcmpi(keys,'keep_heart'));
    if ~isempty(idx)
        params.keep_heart = onoff(varargin{idx+1}, 'keep_heart');
    else
        params.keep_heart = false;
    end

    % get_RR options
    idx = find(strcmpi(keys,'ecg_searchback'));
    if ~isempty(idx)
        params.ecg_searchback = varargin{idx+1};
    end

    % get_RR detection options (see get_RR.m for defaults)
    getRRopts = {'ecg_bandpass' 'ecg_peakthresh' 'ecg_refperiod' 'ecg_polarity' ...
        'ecg_adaptive_pol' 'ppg_bandpass' 'ppg_height_method' 'ppg_polarity'};
    for iOpt = 1:length(getRRopts)
        idx = find(strcmpi(keys,getRRopts{iOpt}));
        if ~isempty(idx)
            params.(getRRopts{iOpt}) = varargin{idx+1};
        end
    end

    % Ectopic beats in the RR cleaning (CLEAN_RR): 'adaptive' (default;
    % Lipponen & Tarvainen, 2019), 'ratio' (fixed 20% rule of v1.6) or 'none'
    idx = find(strcmpi(keys,'ectopic_rule'));
    if ~isempty(idx)
        params.ectopic_rule = lower(char(varargin{idx+1}));
        if ~any(strcmp(params.ectopic_rule, {'adaptive','ratio','none'}))
            error("'ectopic_rule' must be 'adaptive', 'ratio' or 'none'.")
        end
    end

    % Heart components: minimum ICLabel heart probability to remove them
    % ('hep' with 'clean_eeg' and 'heart_removal','ica'; default .75)
    idx = find(strcmpi(keys,'conf_thresh'));
    if ~isempty(idx)
        params.conf_thresh = varargin{idx+1};
    end
    % Heart artifact removal in the HEP cleaning: 'none' (default),
    % 'ica' (ICLabel heart components), 'ica_ecg' (components correlated
    % with the ECG), 'ica_caracas' (CARACAS, SASICA plugin), 'ecg_regression'
    % or 'ssp' (ECG only). Not recommended:
    % each removed little of the artifact or part of a synthetic
    % heartbeat-locked response on 6 ds003838 subjects (see the README)
    idx = find(strcmpi(keys,'heart_removal'));
    if ~isempty(idx)
        params.heart_removal = lower(varargin{idx+1});
        if ~any(strcmp(params.heart_removal, {'none','ica','ica_ecg','ica_caracas','ecg_regression','ssp'}))
            error("'heart_removal' must be 'none', 'ica', 'ica_ecg', 'ica_caracas', 'ecg_regression' or 'ssp'.")
        end
        if any(strcmp(params.heart_removal, {'ica_ecg','ecg_regression','ssp'})) && ~strcmpi(params.heart_signal,'ecg')
            error("'heart_removal','%s' needs an ECG ('heart_signal','ecg').", params.heart_removal)
        end
    else
        params.heart_removal = 'none';
    end
    idx = find(strcmpi(keys,'ecg_corr'));
    if ~isempty(idx), params.ecg_corr = varargin{idx+1}; end
    idx = find(strcmpi(keys,'ssp_n'));
    if ~isempty(idx), params.ssp_n = varargin{idx+1}; end
    % PPG fiducial: detect pulse-wave valleys (default) or peaks (get_RR.m)
    idx = find(strcmpi(keys,'ppg_detect_mode'));
    if ~isempty(idx)
        params.ppg_detect_mode = varargin{idx+1};
        if ~any(strcmpi(params.ppg_detect_mode, {'valleys','peaks'}))
            error("'ppg_detect_mode' must be 'valleys' or 'peaks'.")
        end
    end

    % HEP baseline: no correction (v1.7: the regression option was removed;
    % each heartbeat's baseline is exported in EEG.brainbeats.prestim).
    % 'hep_baseline' is still accepted, with a warning unless 'none'.
    idx = find(strcmpi(keys,'hep_baseline'));
    if ~isempty(idx) && ~strcmpi(varargin{idx+1}, 'none')
        warning('BrainBeats:removed', ['''hep_baseline'' was removed in BrainBeats 1.7: the HEP is not baseline-corrected. ' ...
            'Each heartbeat''s baseline (''hep_baseline_win'') is exported in EEG.brainbeats.prestim.baseline, ' ...
            'for a regression-based correction in the group-level model (Alday, 2019), or use BASELINE_REGRESSION ' ...
            'on epochs that pool conditions.'])
    end
    % Window of each heartbeat's exported baseline voltage (ms)
    idx = find(strcmpi(keys,'hep_baseline_win'));
    if ~isempty(idx)
        params.hep_baseline_win = varargin{idx+1};
        if ~isnumeric(params.hep_baseline_win) || numel(params.hep_baseline_win) ~= 2 || ...
                params.hep_baseline_win(1) >= params.hep_baseline_win(2) || params.hep_baseline_win(2) > 0
            error("'hep_baseline_win' must be [start end] in ms, before the R-peak (e.g. [-150 -50]).")
        end
    end

    % HEP epoch window: [start end] in ms (default [-300 600]) or 'adaptive'
    idx = find(strcmpi(keys,'hep_window'));
    if ~isempty(idx)
        params.hep_window = varargin{idx+1};
        if (ischar(params.hep_window) || isstring(params.hep_window))
            if ~strcmpi(params.hep_window,'adaptive')
                error("'hep_window' must be [start end] in ms or 'adaptive'.")
            end
            params.hep_window = 'adaptive';
        elseif ~isnumeric(params.hep_window) || numel(params.hep_window) ~= 2 || params.hep_window(1) >= 0 || params.hep_window(2) <= 0
            error("'hep_window' must be [start end] in ms (start < 0 < end) or 'adaptive'.")
        end
    end

    % HEP time-frequency measures (HRSP, HRPC) for all channels, and the
    % surrogate heartbeat control (number of surrogates, 0 = none)
    idx = find(strcmpi(keys,'hep_tf_freqs'));
    if ~isempty(idx)
        params.hep_tf_freqs = varargin{idx+1};
        if ~isnumeric(params.hep_tf_freqs) || numel(params.hep_tf_freqs) ~= 2 || params.hep_tf_freqs(1) < 1 || params.hep_tf_freqs(2) <= params.hep_tf_freqs(1)
            error("'hep_tf_freqs' must be [low high] in Hz (e.g. [4 30]).")
        end
    end
    idx = find(strcmpi(keys,'hep_surrogates'));
    if ~isempty(idx)
        params.hep_surrogates = varargin{idx+1};
        if ~isscalar(params.hep_surrogates) || params.hep_surrogates < 0 || params.hep_surrogates ~= round(params.hep_surrogates)
            error("'hep_surrogates' must be the number of surrogates (e.g. 100), or 0 for none.")
        end
    end

    % Surrogate trains: beats at random cardiac phases only (v1.7: 'rigid',
    % the whole train shifted, was removed: still heartbeat-locked, it tests
    % contrasts between conditions, not one condition)
    idx = find(strcmpi(keys,'hep_surrogate_mode'));
    if ~isempty(idx) && ~strcmpi(varargin{idx+1}, 'random')
        warning('BrainBeats:removed', ['''hep_surrogate_mode'' was removed in BrainBeats 1.7: the surrogate beats ' ...
            'are drawn at random cardiac phases. Surrogates made by shifting the whole heartbeat train are still ' ...
            'heartbeat-locked, so they test contrasts between conditions, not one condition.'])
    end
    idx = find(strcmpi(keys,'hep_surr_stat'));
    if ~isempty(idx)
        params.hep_surr_stat = lower(varargin{idx+1});
        if ~any(strcmp(params.hep_surr_stat, {'tfce','tmax','fdr'}))
            error("'hep_surr_stat' must be 'tfce', 'tmax' or 'fdr'.")
        end
    end
    % HRSP plotted: power with the HEP subtracted from each heartbeat
    % ('induced', default) or all power ('total'); both are stored
    idx = find(strcmpi(keys,'hep_tf_power'));
    if ~isempty(idx)
        params.hep_tf_power = lower(varargin{idx+1});
        if ~any(strcmp(params.hep_tf_power, {'induced','total'}))
            error("'hep_tf_power' must be 'induced' or 'total'.")
        end
    end
    % QRS map of the cardiac field artifact and its stability (default on)
    idx = find(strcmpi(keys,'hep_cfa'));
    if ~isempty(idx)
        params.hep_cfa = onoff(varargin{idx+1}, 'hep_cfa');
    end

    % HEP measures on the scalp channels, the independent components, or
    % both; ROI of the channel plots
    idx = find(strcmpi(keys,'hep_level'));
    if ~isempty(idx)
        params.hep_level = lower(varargin{idx+1});
        if ~any(strcmp(params.hep_level, {'channels','ics','both'}))
            error("'hep_level' must be 'channels', 'ics' or 'both'.")
        end
    end
    idx = find(strcmpi(keys,'hep_roi'));
    if ~isempty(idx)
        params.hep_roi = cellstr(varargin{idx+1});
    end

    % PPG pulse arrival time, to shift PPG beats back to the heartbeats for
    % HEP: 'auto' (default: from an ECG channel of the file if there is one,
    % else an error that suggests the estimate from the EEG cardiac field
    % artifact), the label of an ECG channel, 'eeg', a delay in ms, or 'off'
    idx = find(strcmpi(keys,'ppg_transit'));
    if ~isempty(idx)
        params.ppg_transit = varargin{idx+1};
        if iscell(params.ppg_transit), params.ppg_transit = params.ppg_transit{1}; end
        if ~strcmpi(params.heart_signal,'ppg')
            error("'ppg_transit' only applies to 'heart_signal','ppg'.")
        end
    elseif strcmpi(params.heart_signal,'ppg') && strcmp(params.analysis,'hep')
        params.ppg_transit = 'auto';
    end
end

%% HRV features

if params.heart && strcmp(params.analysis,'features')
    idx = find(strcmpi(keys,'hrv_features'));
    if ~isempty(idx) && ((ischar(varargin{idx+1}) && strcmpi(varargin{idx+1},'off')) || isequal(varargin{idx+1},false))
        params.hrv_features = false;
        params.hrv_time = false;
        params.hrv_frequency = false;
        params.hrv_nonlinear = false;
    elseif ~isempty(idx)
        if ~isfield(params,'hrv_features')
            params.hrv_features = true;
        end
        hrv_features = varargin{idx+1};
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
        % All domains by default
        disp('HRV features not defined. Setting all domains by default (time, frequency, and nonlinear)')
        params.hrv_features = true;
        params.hrv_time = true;
        params.hrv_frequency = true;
        params.hrv_nonlinear = true;
    end
else
    params.hrv_features = false;
    params.hrv_time = false;
    params.hrv_frequency = false;
    params.hrv_nonlinear = false;
end

if params.hrv_frequency
    idx = find(strcmpi(keys,'hrv_norm'));
    if ~isempty(idx)
        params.hrv_norm = varargin{idx+1};
    end
    idx = find(strcmpi(keys,'hrv_spec'));
    if ~isempty(idx)
        params.hrv_spec = varargin{idx+1};
    end
    idx = find(strcmpi(keys,'hrv_overlap'));
    if ~isempty(idx)
        params.hrv_overlap = varargin{idx+1};
    end
end
idx = find(strcmpi(keys,'prsa_thresh'));
if params.hrv_nonlinear && ~isempty(idx)
    params.prsa_thresh = varargin{idx+1};
end

%% EEG parameters

% Any EEG operations? ('eeg','off' turns them all off)
idx = find(strcmpi(keys,'eeg'));
if ~isempty(idx)
    params.eeg = varargin{idx+1};
    if strcmpi(params.eeg,'off')
        params.eeg = false;
        fprintf('Turning OFF all EEG operations. \n')
    end
else
    params.eeg = true;
end

% Preprocess EEG
if params.eeg
    idx = find(strcmpi(keys,'clean_eeg'));
    if ~isempty(idx)
        params.clean_eeg = onoff(varargin{idx+1}, 'clean_eeg');
    else
        warning("'clean_eeg' input not defined. Using default: NO preprocessing (i.e., assuming you have already preprocessed your EEG data. If you wish to preprocess your EEG data with BrainBeats, set 'clean_eeg' to true.")
        params.clean_eeg = false;
    end

    % EEG preprocessing
    idx = find(strcmpi(keys,'highpass'));
    if ~isempty(idx)
        params.highpass = varargin{idx+1};
    end
    idx = find(strcmpi(keys,'lowpass'));
    if ~isempty(idx)
        params.lowpass = varargin{idx+1};
    end
    idx = find(strcmpi(keys,'flatline'));
    if ~isempty(idx)
        params.flatline = varargin{idx+1};
    end
    idx = find(strcmpi(keys,'corrThresh'));
    if ~isempty(idx)
        params.corrThresh = varargin{idx+1};
    end
    idx = find(strcmpi(keys,'maxBad'));
    if ~isempty(idx)
        params.maxBad = varargin{idx+1};
    end
    idx = find(strcmpi(keys,'asr_cutoff'));
    if ~isempty(idx)
        params.asr_cutoff = varargin{idx+1};
    end
    if any(strcmpi(keys,'asr_mem'))
        warning("'asr_mem' is ignored since v1.7: ASR runs in its default chunks, as its output depends on the chunk size.")
    end
    idx = find(strcmpi(keys,'ref'));
    if ~isempty(idx)
        params.ref = lower(varargin{idx+1});
        if ~any(strcmp(params.ref, {'average','infinity','csd','off'}))
            error("'ref' must be 'average', 'infinity', 'csd' (surface Laplacian) or 'off'.")
        end
    end
    % Artifact removal: 'asr_ica' (default: ASR or bad epochs, then ICA +
    % ICLabel) or 'gedai' (GEDAI plugin, installed if missing)
    idx = find(strcmpi(keys,'clean_method'));
    if ~isempty(idx)
        params.clean_method = lower(varargin{idx+1});
        if ~any(strcmp(params.clean_method, {'asr_ica','gedai'}))
            error("'clean_method' must be 'asr_ica' or 'gedai'.")
        end
    end
    idx = find(strcmpi(keys,'linenoise'));
    if ~isempty(idx)
        params.linenoise = varargin{idx+1};
    end
    idx = find(strcmpi(keys,'filttype'));
    if ~isempty(idx)
        params.filttype = varargin{idx+1};
    end
    idx = find(strcmpi(keys,'detectMethod'));
    if ~isempty(idx)
        params.detectMethod = varargin{idx+1};
    end
    idx = find(strcmpi(keys,'icamethod') | strcmpi(keys,'ica_method'));
    if ~isempty(idx)
        params.icamethod = varargin{idx(1)+1};
        if ~isscalar(params.icamethod) || ~any(params.icamethod == [1 2 3])
            error("'icamethod' must be 1 (Picard), 2 (Infomax), or 3 (replicable Infomax).")
        end
    end
    % Data the HEP ICA is fitted on: a continuous 1-Hz high-passed copy,
    % with ASR at this cutoff ('ica_asr', default 80, 0 = off), this
    % low-pass in Hz ('ica_lowpass', default 0 = none) and sampling rate
    % ('ica_srate', default: that of the data). See CLEAN_EEG.
    for k = {'ica_asr' 'ica_lowpass' 'ica_srate'}
        idx = find(strcmpi(keys,k{1}));
        if ~isempty(idx)
            v = varargin{idx+1};
            if ~isnumeric(v) || ~isscalar(v) || v < 0 || (v == 0 && strcmp(k{1},'ica_srate'))
                error("'%s' must be a positive number ('ica_asr' and 'ica_lowpass' can be 0 for none).", k{1})
            end
            params.(k{1}) = v;
        end
    end
    if isfield(params,'ica_lowpass') && params.ica_lowpass > 0 && params.ica_lowpass <= 2
        error("'ica_lowpass' must be 0 (none) or above the 1-Hz high-pass of the ICA copy (e.g. 30).")
    end
    if isfield(params,'ica_srate') && (params.ica_srate < 50 || ...
            (isfield(params,'ica_lowpass') && params.ica_lowpass > 0 && params.ica_srate < 2*params.ica_lowpass))
        error("'ica_srate' must be at least 50 Hz and twice 'ica_lowpass'.")
    end
    idx = find(strcmpi(keys,'ica_epochs'));   % epoch that copy like the HEP data (bad epochs removed from it)
    if ~isempty(idx)
        params.ica_epochs = onoff(varargin{idx+1}, 'ica_epochs');
    end

    % Cleaning options have no effect without 'clean_eeg'
    if ~params.clean_eeg
        ignored = {'lowpass' 'flatline' 'corrThresh' 'maxBad' 'asr_cutoff' 'asr_mem' ...
            'ref' 'clean_method' 'linenoise' 'detectMethod' 'heart_removal' 'ecg_corr' 'ssp_n' ...
            'ica_asr' 'ica_lowpass' 'ica_srate' 'ica_epochs'};
        if ~isfield(params,'hep_level') || strcmp(params.hep_level,'channels')
            ignored = [ignored {'highpass' 'filttype'}];    % else used by the ICA of the IC measures
        end
        ignored = given(ismember(lower(given), lower(ignored)));
        if ~isempty(ignored)
            warning("Input(s) ignored because 'clean_eeg' is false: %s.", strjoin(ignored, ', '))
        end
    end
    % HEP options have no effect on the other analyses
    if ~strcmp(params.analysis,'hep')
        ignored = given(startsWith(lower(given), {'ica_' 'hep_'}) & ~strcmpi(given, 'ica_method') | ...
            ismember(lower(given), {'keep_heart' 'ppg_transit'}));   % 'ica_method': alias of 'icamethod', used by any cleaning
        if ~isempty(ignored)
            warning("Input(s) ignored for '%s' (HEP only): %s.", params.analysis, strjoin(ignored, ', '))
        end
    end
end

%% EEG features

% Turn OFF all EEG features
if ~params.eeg
    params.eeg_features = false;
    params.eeg_time = false;
    params.eeg_frequency = false;
    params.eeg_nonlinear = false;
    params.clean_eeg = false;
end

if params.eeg
    if strcmpi(params.analysis,'features')

        params.eeg_features = true;

        % Which features to compute ('off' for none)
        idx = find(strcmpi(keys,'eeg_features'));
        if ~isempty(idx) && ((ischar(varargin{idx+1}) && strcmpi(varargin{idx+1},'off')) || isequal(varargin{idx+1},false))
            params.eeg_features = false;
            params.eeg_time = false;
            params.eeg_frequency = false;
            params.eeg_nonlinear = false;
        elseif ~isempty(idx)
            eeg_features = varargin{idx+1};
            if sum(contains(eeg_features,{'time'})) > 0
                params.eeg_time = true;
            else
                params.eeg_time = false;
            end
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
            disp('EEG features not specified. Setting all EEG features (time, frequency, nonlinear) to ON. ')
            params.eeg_features = true;
            params.eeg_time = true;
            params.eeg_frequency = true;
            params.eeg_nonlinear = true;
        end

        % EEG frequency-domain parameters
        if params.eeg_frequency
            idx = find(strcmpi(keys,'eeg_frange'));
            if ~isempty(idx)
                params.eeg_frange = varargin{idx+1};
            end
            idx = find(strcmpi(keys,'eeg_wintype'));
            if ~isempty(idx)
                params.eeg_wintype = varargin{idx+1};
            end
            idx = find(strcmpi(keys,'eeg_winoverlap'));
            if ~isempty(idx)
                params.eeg_winoverlap = varargin{idx+1};
            end
            idx = find(strcmpi(keys,'eeg_winlen'));
            if ~isempty(idx)
                params.eeg_winlen = varargin{idx+1};
            end
            idx = find(strcmpi(keys,'eeg_freqbounds'));
            if ~isempty(idx)
                params.eeg_freqbounds = varargin{idx+1};
            end
            idx = find(strcmpi(keys,'eeg_norm'));
            if ~isempty(idx)
                params.eeg_norm = varargin{idx+1};
            end
            idx = find(strcmpi(keys,'asy_norm'));
            if ~isempty(idx)
                params.asy_norm = varargin{idx+1};
            end
        end
    end
end % if params.eeg

%% Visualization and saving

% Visualize preprocessings
idx = find(strcmpi(keys,'vis_cleaning'));
if ~isempty(idx)
    params.vis_cleaning = onoff(varargin{idx+1}, 'vis_cleaning');
    if params.vis_cleaning
        fprintf('Visualization of data cleaning set to ON. \n')
    else
        fprintf('Visualization of data cleaning set to OFF. \n')
    end
else
    params.vis_cleaning = true;
    disp("Visualization of data cleaning not defined. Set to ON by default. If you wish to turn it OFF, set input 'vis_cleaning' to false");
end

% Visualize outputs
idx = find(strcmpi(keys,'vis_outputs'));
if ~isempty(idx)
    params.vis_outputs = onoff(varargin{idx+1}, 'vis_outputs');
    if params.vis_outputs
        fprintf('Visualization of outputs set to ON. \n')
    else
        fprintf('Visualization of outputs set to OFF. \n')
    end
else
    params.vis_outputs = true;
    disp("Visualization of outputs not defined. Set to ON by default. If you wish to turn it OFF, set input 'vis_outputs' to false");
end

% Save outputs
idx = find(strcmpi(keys,'save'));
if ~isempty(idx)
    params.save = onoff(varargin{idx+1}, 'save');
    if params.save
        fprintf('Saving outputs set to ON. \n')
    else
        fprintf('Saving outputs set to OFF. \n')
    end
else
    params.save = true; % ON by default
    disp("Saving outputs not defined. Set to ON by default. If you wish to turn it OFF, set input 'save' to false");
end

%% Parallel & GPU computing

% Parallel computing
idx = find(strcmpi(keys,'parpool'));
if ~isempty(idx)
    params.parpool = onoff(varargin{idx+1}, 'parpool');
else
    params.parpool = false;
    fprintf('Parallel computing not defined: set to OFF by default. \n')
end


% GPU computing
idx = find(strcmpi(keys,'gpu'));
if ~isempty(idx)
    params.gpu = onoff(varargin{idx+1}, 'gpu');
    if params.gpu, fprintf('GPU computing: set to ON. \n'); end
else
    params.gpu = false;
    fprintf('GPU computing not defined: set to OFF by default. \n')
end

idx = find(strcmpi(keys,'gong'));
if ~isempty(idx)
    params.gong = onoff(varargin{idx+1}, 'gong');
else
    params.gong = true;
end

function v = onoff(v, name)
% Logical scalar from true/false, 1/0, 'on'/'off', 'true'/'false', 'yes'/'no'
if ischar(v) || isstring(v)
    s = lower(char(v));
    if any(strcmp(s, {'on' 'true' 'yes'})), v = true; return, end
    if any(strcmp(s, {'off' 'false' 'no'})), v = false; return, end
end
if isscalar(v) && (islogical(v) || (isnumeric(v) && ~isnan(v)))
    v = logical(v);
else
    error('BrainBeats:input', "'%s' must be true or false (or 'on'/'off').", name)
end

function k = known_keys()
% Inputs of BRAINBEATS_PROCESS
k = {'heart_signal' 'heart_channels' 'beat_latencies' 'analysis' 'keep_heart' ...
    'ecg_searchback' 'ecg_bandpass' 'ecg_peakthresh' 'ecg_refperiod' 'ecg_polarity' ...
    'ecg_adaptive_pol' 'ppg_bandpass' 'ppg_height_method' 'ppg_polarity' 'ppg_detect_mode' 'ppg_transit' 'ectopic_rule' ...
    'conf_thresh' 'heart_removal' 'ecg_corr' 'ssp_n' 'hep_baseline' 'hep_baseline_win' 'hep_window' ...
    'hep_tf_freqs' 'hep_surrogates' 'hep_surrogate_mode' 'hep_surr_stat' 'hep_tf_power' 'hep_cfa' 'hep_level' 'hep_roi' ...
    'hrv_features' 'hrv_norm' 'hrv_spec' 'hrv_overlap' 'prsa_thresh' 'eeg' 'clean_eeg' 'highpass' 'lowpass' ...
    'flatline' 'corrThresh' 'maxBad' 'asr_cutoff' 'asr_mem' 'ref' 'clean_method' 'linenoise' ...
    'filttype' 'detectMethod' 'icamethod' 'ica_method' 'ica_asr' 'ica_lowpass' 'ica_srate' 'ica_epochs' ...
    'eeg_features' 'eeg_frange' ...
    'eeg_wintype' 'eeg_winoverlap' 'eeg_winlen' 'eeg_freqbounds' 'eeg_norm' 'asy_norm' ...
    'vis_cleaning' 'vis_outputs' 'save' 'parpool' 'gpu' 'gong'};
k = [k legacy_keys()];

function k = legacy_keys()
% Inputs of earlier versions, accepted and ignored
k = {'ppg_learnperiod' 'ppg_buffer' 'ppg_learnthresh' 'ppg_eyeclosing' 'ppg_expctperiod' ...
    'ppg_slopewindow' 'rr_physlimlow' 'rr_physlimhigh' 'rr_gaplim' 'rr_changelim' 'rr_correct' ...
    'clean_heart' 'eeg_interp' 'boost'};
