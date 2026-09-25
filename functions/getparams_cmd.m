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

%% General parameters

% Extract user parameters
idx = find(strcmpi(varargin,'heart_signal'));
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
    idx = find(strcmpi(varargin,'heart_channels'));
    if ~isempty(idx)
        params.heart_channels = varargin{idx+1};
        fprintf('Number of heart channels selected: %g \n', length(params.heart_channels));
    elseif ~strcmpi(params.heart_signal, 'rr')
        error("Heart channels not defined. Please define 'heart_channels'. See help for an example")
    else
        params.heart_channels = {};
    end

    % Pre-detected beat latencies (for 'rr' mode)
    idx = find(strcmpi(varargin,'beat_latencies'));
    if ~isempty(idx)
        params.beat_latencies = varargin{idx+1};
        fprintf('Pre-detected beat latencies provided: %g beats \n', length(params.beat_latencies));
    elseif strcmpi(params.heart_signal, 'rr')
        error("'beat_latencies' required when heart_signal is 'rr'. Provide a vector of beat times in seconds.")
    end
end

% Analysis to do
idx = find(strcmpi(varargin,'analysis'));
if ~isempty(idx)
    params.analysis = lower(varargin{idx+1});
    fprintf('Selected mode to run: %s \n', params.analysis);
else
    error("Analysis type not defined. Must be set to either 'hep', 'features', 'rm_heart', or 'coherence'. ")
end
if ~any(strcmp(params.analysis, {'hep' 'features' 'rm_heart' 'coherence'}))
    error("Analysis '%s' not recognized. Must be 'hep', 'features', 'rm_heart', or 'coherence'.", params.analysis)
end
if ~params.heart && any(strcmp(params.analysis, {'hep' 'rm_heart' 'coherence'}))
    error("The '%s' analysis needs a heart signal: 'heart_signal' cannot be 'off'.", params.analysis)
end
if strcmp(params.heart_signal,'rr') && any(strcmp(params.analysis, {'rm_heart' 'coherence'}))
    error("The '%s' analysis needs the ECG/PPG signal itself: 'heart_signal' cannot be 'rr'.", params.analysis)
end

%% Cardiovascular signal preprocessing

if params.heart
    idx = find(strcmpi(varargin,'clean_heart'));
    if ~isempty(idx)
        params.clean_heart = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'keep_heart'));
    if ~isempty(idx)
        params.keep_heart = varargin{idx+1};
    else
        params.keep_heart = false;
    end

    % Legacy PPG detector options (kept for compatibility, not used by get_RR)
    idx = find(strcmpi(varargin,'ppg_learnperiod'));
    if ~isempty(idx)
        params.ppg_learnperiod = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'ppg_buffer'));
    if ~isempty(idx)
        params.ppg_buffer = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'ppg_learnthresh'));
    if ~isempty(idx)
        params.ppg_learnthresh = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'ppg_eyeclosing'));
    if ~isempty(idx)
        params.ppg_eyeclosing = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'ppg_expctperiod'));
    if ~isempty(idx)
        params.ppg_expctperiod = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'ppg_slopewindow'));
    if ~isempty(idx)
        params.ppg_slopewindow = varargin{idx+1};
    end

    % get_RR options
    idx = find(strcmpi(varargin,'ecg_searchback'));
    if ~isempty(idx)
        params.ecg_searchback = varargin{idx+1};
    end

    % get_RR detection options (see get_RR.m for defaults)
    getRRopts = {'ecg_bandpass' 'ecg_peakthresh' 'ecg_refperiod' 'ecg_polarity' ...
        'ecg_adaptive_pol' 'ppg_bandpass' 'ppg_height_method'};
    for iOpt = 1:length(getRRopts)
        idx = find(strcmpi(varargin,getRRopts{iOpt}));
        if ~isempty(idx)
            params.(getRRopts{iOpt}) = varargin{idx+1};
        end
    end

    % params for method 3: removing heart artifacts from EEG
    idx = find(strcmpi(varargin,'conf_thresh'));
    if ~isempty(idx)
        params.conf_thresh = varargin{idx+1};
    end
    if any(strcmpi(varargin,'boost'))
        warning("The 'boost' option was removed (it did not improve the detection of heart components). Ignoring it.")
    end


    % Legacy RR-correction options (kept for compatibility, not used by clean_rr)
    idx = find(strcmpi(varargin,'rr_physlimlow'));
    if ~isempty(idx)
        params.rr_physlimlow = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'rr_physlimhigh'));
    if ~isempty(idx)
        params.rr_physlimhigh = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'rr_gaplim'));
    if ~isempty(idx)
        params.rr_gaplim = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'rr_changelim'));
    if ~isempty(idx)
        params.rr_changelim = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'rr_correct'));
    if ~isempty(idx)
        params.rr_correct = varargin{idx+1};
    end

    % PPG fiducial: detect pulse-wave valleys (default) or peaks (get_RR.m)
    idx = find(strcmpi(varargin,'ppg_detect_mode'));
    if ~isempty(idx)
        params.ppg_detect_mode = varargin{idx+1};
        if ~any(strcmpi(params.ppg_detect_mode, {'valleys','peaks'}))
            error("'ppg_detect_mode' must be 'valleys' or 'peaks'.")
        end
    end

    % HEP baseline: 'none' (default) or 'regression' (Alday, 2019)
    idx = find(strcmpi(varargin,'hep_baseline'));
    if ~isempty(idx)
        params.hep_baseline = lower(varargin{idx+1});
        if ~any(strcmp(params.hep_baseline, {'none','regression'}))
            error("'hep_baseline' must be 'none' or 'regression'.")
        end
    end
    idx = find(strcmpi(varargin,'hep_baseline_win'));
    if ~isempty(idx)
        params.hep_baseline_win = varargin{idx+1};
    end

    % HEP epoch window: [start end] in ms (default [-300 600]) or 'adaptive'
    idx = find(strcmpi(varargin,'hep_window'));
    if ~isempty(idx)
        params.hep_window = varargin{idx+1};
        if (ischar(params.hep_window) || isstring(params.hep_window))
            if ~strcmpi(params.hep_window,'adaptive')
                error("'hep_window' must be [start end] in ms or 'adaptive'.")
            end
            params.hep_window = 'adaptive';
        elseif ~isnumeric(params.hep_window) || numel(params.hep_window) ~= 2 || params.hep_window(1) >= 0
            error("'hep_window' must be [start end] in ms (start < 0) or 'adaptive'.")
        end
    end

    % PPG pulse arrival time, to shift PPG beats back to the heartbeat for
    % HEP: a delay in ms, or the label of an ECG channel to estimate it from
    idx = find(strcmpi(varargin,'ppg_transit'));
    if ~isempty(idx)
        params.ppg_transit = varargin{idx+1};
        if iscell(params.ppg_transit), params.ppg_transit = params.ppg_transit{1}; end
        if ~strcmpi(params.heart_signal,'ppg')
            error("'ppg_transit' only applies to 'heart_signal','ppg'.")
        end
    end
end

%% HRV features

if params.heart && strcmp(params.analysis,'features')
    idx = find(strcmpi(varargin,'hrv_features'));
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
    idx = find(strcmpi(varargin,'hrv_norm'));
    if ~isempty(idx)
        params.hrv_norm = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'hrv_spec'));
    if ~isempty(idx)
        params.hrv_spec = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'hrv_overlap'));
    if ~isempty(idx)
        params.hrv_overlap = varargin{idx+1};
    end
end

%% EEG parameters

% Any EEG operations? ('eeg','off' turns them all off)
idx = find(strcmpi(varargin,'eeg'));
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
    idx = find(strcmpi(varargin,'clean_eeg'));
    if ~isempty(idx)
        params.clean_eeg = varargin{idx+1};
    else
        warning("'clean_eeg' input not defined. Using default: NO preprocessing (i.e., assuming you have already preprocessed your EEG data. If you wish to preprocess your EEG data with BrainBeats, set 'clean_eeg' to true.")
        params.clean_eeg = false;
    end

    % EEG preprocessing
    idx = find(strcmpi(varargin,'highpass'));
    if ~isempty(idx)
        params.highpass = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'lowpass'));
    if ~isempty(idx)
        params.lowpass = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'flatline'));
    if ~isempty(idx)
        params.flatline = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'corrThresh'));
    if ~isempty(idx)
        params.corrThresh = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'maxBad'));
    if ~isempty(idx)
        params.maxBad = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'eeg_interp'));
    if ~isempty(idx)
        params.eeg_interp = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'asr_cutoff'));
    if ~isempty(idx)
        params.asr_cutoff = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'asr_mem'));
    if ~isempty(idx)
        params.asr_mem = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'ref'));
    if ~isempty(idx)
        params.ref = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'linenoise'));
    if ~isempty(idx)
        params.linenoise = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'filttype'));
    if ~isempty(idx)
        params.filttype = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'detectMethod'));
    if ~isempty(idx)
        params.detectMethod = varargin{idx+1};
    end
    idx = find(strcmpi(varargin,'icamethod') | strcmpi(varargin,'ica_method'));
    if ~isempty(idx)
        params.icamethod = varargin{idx(1)+1};
        if ~isscalar(params.icamethod) || ~any(params.icamethod == [1 2 3])
            error("'icamethod' must be 1 (Picard), 2 (Infomax), or 3 (replicable Infomax).")
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
        idx = find(strcmpi(varargin,'eeg_features'));
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
            idx = find(strcmpi(varargin,'eeg_frange'));
            if ~isempty(idx)
                params.eeg_frange = varargin{idx+1};
            end
            idx = find(strcmpi(varargin,'eeg_wintype'));
            if ~isempty(idx)
                params.eeg_wintype = varargin{idx+1};
            end
            idx = find(strcmpi(varargin,'eeg_winoverlap'));
            if ~isempty(idx)
                params.eeg_winoverlap = varargin{idx+1};
            end
            idx = find(strcmpi(varargin,'eeg_winlen'));
            if ~isempty(idx)
                params.eeg_winlen = varargin{idx+1};
            end
            idx = find(strcmpi(varargin,'eeg_freqbounds'));
            if ~isempty(idx)
                params.eeg_freqbounds = varargin{idx+1};
            end
            idx = find(strcmpi(varargin,'eeg_norm'));
            if ~isempty(idx)
                params.eeg_norm = varargin{idx+1};
            end
            idx = find(strcmpi(varargin,'asy_norm'));
            if ~isempty(idx)
                params.asy_norm = varargin{idx+1};
            end
        end
    end
end % if params.eeg

%% Visualization and saving

% Visualize preprocessings
idx = find(strcmpi(varargin,'vis_cleaning'));
if ~isempty(idx)
    params.vis_cleaning = varargin{idx+1};
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
idx = find(strcmpi(varargin,'vis_outputs'));
if ~isempty(idx)
    params.vis_outputs = varargin{idx+1};
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
idx = find(strcmpi(varargin,'save'));
if ~isempty(idx)
    params.save = varargin{idx+1};
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
idx = find(strcmpi(varargin,'parpool'));
if ~isempty(idx)
    val = varargin{idx+1};
    if ischar(val) || isstring(val)
        val = strcmpi(val,'on');   % 'on' / 'off'
    end
    params.parpool = logical(val);
else
    params.parpool = false;
    fprintf('Parallel computing not defined: set to OFF by default. \n')
end


% GPU computing
idx = find(strcmpi(varargin,'gpu'));
if ~isempty(idx)
    params.gpu = logical(varargin{idx+1});
    if params.gpu, fprintf('GPU computing: set to ON. \n'); end
else
    params.gpu = false;
    fprintf('GPU computing not defined: set to OFF by default. \n')
end

idx = find(strcmpi(varargin,'gong'));
if ~isempty(idx)
    params.gong = varargin{idx+1};
else
    params.gong = true;
end
