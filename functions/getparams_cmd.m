%% Extract all the inputs and parameters specified by user via command line
%
% Cedric Cannard, April 2023

function params = getparams_cmd(varargin)

%% General parameters

% Extract user parameters
idx = find(strcmpi(varargin,'heart_signal'));
if ~isempty(idx)
    params.heart = true;
    params.heart_signal = lower(varargin{idx+1});
    if ~contains(params.heart_signal, {'ecg' 'ppg' 'off'})
        error("Heart signal not recognized. Should be 'ecg' or 'ppg' or 'off'.")
    end
    if contains(params.heart_signal, 'off')
        params.heart = false;
        params.heart_channels = {};
    end
else
    error("Heart signal type not defined. Please define 'heart_signal' as 'ecg', 'ppg', or 'off'. Type help in the command window for an example")
end

% Detect if there are 1 or multiple heart channel(s)
if params.heart
    idx = find(strcmpi(varargin,'heart_channels'));
    if ~isempty(idx)
        params.heart_channels = varargin{idx+1};
        fprintf('Number of heart channels selected: %g \n', length(params.heart_channels));
    else
        error("Heart channels not defined. Please define 'heart_channels'. See help for an example")
    end
end

% Analysis to do
idx = find(strcmpi(varargin,'analysis'));
if ~isempty(idx)
    params.analysis = varargin{idx+1};
    fprintf('Selected mode to run: %s \n', params.analysis);
else
    error("Analysis type not defined. Must be set to either 'hep', 'features', or 'rm_heart'. ")
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

    % get_RR on PPG signal
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

    % get_RR on ECG signal
    idx = find(strcmpi(varargin,'ecg_searchback'));
    if ~isempty(idx)
        params.ecg_searchback = varargin{idx+1};
    end

    % params for method 3:removing heart artifacts from EEG
    if isfield(params, 'conf_thresh') && ~isempty(params.conf_thresh)
        params.conf_thresh = params.conf_thresh;
    end
    if isfield(params, 'boost') && ~isempty(params.boost)
        params.boost = logical(params.boost);
    end


    % RR artifacts
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
end

%% HRV features

if params.heart && strcmp(params.analysis,'features')
    idx = find(strcmpi(varargin,'hrv_features'));
    if ~isempty(idx)
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
        % run all domains by default if not set
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

% Any EEG operations?
idx = find(strcmpi(varargin,'eeg'));
if ~isempty(idx)
    params.eeg = varargin{idx+1};  % set 'eeg' to false to turn off all EEG operations
    if strcmpi(params.eeg,'off')
        params.eeg = false;
        fprintf('Turning OFF all EEG operations. \n')
    end
else
    params.eeg = true;  % for any input containing 'eeg' 
end

% Preprocess EEG
if params.eeg
    idx = find(strcmpi(varargin,'clean_eeg'));
    if ~isempty(idx)
        params.clean_eeg = varargin{idx+1};
        % if ~islogical(params.clean_eeg), error("The 'clean_eeg' input should be logical (true or false)."); end
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

        params.eeg_features = true;  % EEG features turned ON

        % which features to compute
        idx = find(strcmpi(varargin,'eeg_features'));
        if ~isempty(idx)
    
    
            % which features to compute
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
    % if ~islogical(params.vis_cleaning), error("The 'vis_cleaning' input should be logical (true or false)."); end
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
    % if ~islogical(params.vis_outputs), error("The 'vis_outputs' input should be logical (true or false)."); end
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
    % if ~islogical(params.save), error("The 'save' input should be logical (true or false)."); end
    if params.save
        fprintf('Saving outputs set to ON. \n')
    else
        fprintf('Saving outputs set to OFF. \n')
    end
else
    params.save = true; % save output by default
    disp("Saving outputs not defined. Set to ON by default. If you wish to turn it OFF, set input 'save' to false");
end

%% Parallel & GPU computing

% Parallel computing
idx = find(strcmpi(varargin,'parpool'));
if ~isempty(idx)
    params.parpool = varargin{idx+1};
    % if ~islogical(params.parpool), error("The 'parpool' input should be logical (true or false)."); end
    if strcmpi(params.parpool,'on') || params.parpool
        params.parpool = true;
    else
        params.parpool = false;
    end
else
    params.parpool = false;
    fprintf('Parallel computing not defined: set to OFF by default. \n')
end


% GPU computing
idx = find(strcmpi(varargin,'gpu'));
if ~isempty(idx)
    params.gpu = varargin{idx+1};
    fprintf('GPU computing: set to ON. \n')
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
