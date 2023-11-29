%% Extract all the inputs and parameters specified by user via command line
% 
% Cedric Cannard, April 2023

function params = getparams_command(varargin)

% Extract user parameters
idx = find(strcmpi(varargin,'heart_signal'));
if ~isempty(idx)
    params.heart_signal = lower(varargin{idx+1});
    if ~contains(params.heart_signal, {'ecg' 'ppg' 'off'})
        error("Heart signal not recognized. Should be 'ecg' or 'ppg' or 'off'.")
    end
    if contains(params.heart_signal, 'off')
        params.heart = false;
    end
else
    error("Heart signal type not defined. Please define 'heart_signal' as 'ecg', 'ppg', or 'off'. Type help in the command window for an example")
end

% Detect if there are 1 or multiple heart channel(s)
idx = find(strcmpi(varargin,'heart_channels'));
if ~isempty(idx)
    params.heart_channels = varargin{idx+1};
    fprintf('Number of heart channels selected: %g \n', length(params.heart_channels));
else
    error("Heart channels not defined. Please define 'heart_channels'. See help for an example")
end

% Analysis to do
idx = find(strcmpi(varargin,'analysis'));
if ~isempty(idx)
    params.analysis = varargin{idx+1};
    fprintf('Selected mode to run: %s \n', params.analysis);
else
    error("Analysis type not defined. Must be set to either 'hep', 'features', or 'rm_heart'. ")
end

% RR artifact correction method
idx = find(strcmpi(varargin,'rr_correct'));
if ~isempty(idx)
    params.rr_correct = varargin{idx+1};
    fprintf('RR artifact method selected by user: %s \n', params.rr_correct)
elseif isempty(idx) && contains(params.analysis, {'features' 'hep'})
    fprintf("RR artifact correction method not defined. Selecting default method: shape-preserving piecewise cubic interpolation ('pchip'). \n")
    params.rr_correct = 'pchip';  % pchip or spline should be default
end

% Clean EEG
idx = find(strcmpi(varargin,'clean_eeg'));
if ~isempty(idx)
    params.clean_eeg = varargin{idx+1};
    if ~islogical(params.clean_eeg), error("The 'clean_eeg' input should be logical (true or false)."); end
else
    warning("'clean_eeg' input not defined. Using default: NO preprocessing (i.e., assuming you have already preprocessed your EEG data. If you wish to preprocess your EEG data with BrainBeats, set 'clean_eeg' to true.")
    params.clean_eeg = false;
end

% EEG features
if strcmp(params.analysis,'features')
    idx = find(strcmpi(varargin,'eeg_features'));
    if ~isempty(idx)
        params.eeg = true;
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
        params.eeg = false;
        fprintf('You chose NOT to extract EEG features on these data. \n')
    end
end

% HRV features
if strcmp(params.analysis,'features')
    idx = find(strcmpi(varargin,'hrv_features'));
    if ~isempty(idx)
        params.hrv = true;
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
        params.hrv = true;
    end
else
    params.hrv = false;
    fprintf('You chose NOT to extract HRV features on these data. \n')
end

% Normalize frequency features
idx = find(strcmpi(varargin,'norm'));
if ~isempty(idx)
    params.norm = varargin{idx+1};
    if ~islogical(params.norm), error("The 'norm' input should be logical (true or false)."); end
    if params.norm
        fprintf('Normalization of frequency-domain outputs set to ON. \n')
    else
        fprintf('Normalization of frequency-domain outputs set to OFF. \n')
    end
else
    params.norm = true;
    fprintf('Normalization of frequency-domain outputs not defined: set to ON by default. \n')
end


% hrv_spec
idx = find(strcmpi(varargin,'hrv_spec'));
if ~isempty(idx)
    params.hrv_spec = varargin{idx+1};
    fprintf('Method to estimate HRV power spectrum set to: %s \n', params.hrv_spec)
else
    params.hrv_spec = 'LombScargle';  % 'LombScargle' (default), 'pwelch', 'fft'
    fprintf("Method to estimate HRV power spectrum not defined: set to 'LombScargle' by default. \n")
end

% hrv_overlap
idx = find(strcmpi(varargin,'hrv_overlap'));
if ~isempty(idx)
    params.hrv_overlap = varargin{idx+1};
    fprintf('Overlap for HRV power spectrum set to: %g%%. \n', params.hrv_overlap*100)
else
    params.hrv_overlap = 0.25;
    fprintf('Overlap for HRV power spectrum not defined and set to default: 25%% overlap \n')
end

% Visualize preprocessings
idx = find(strcmpi(varargin,'vis_cleaning'));
if ~isempty(idx)
    params.vis_cleaning = varargin{idx+1};
    if ~islogical(params.vis_cleaning), error("The 'vis_cleaning' input should be logical (true or false)."); end
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
    if ~islogical(params.vis_outputs), error("The 'vis_outputs' input should be logical (true or false)."); end
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
    if ~islogical(params.save), error("The 'save' input should be logical (true or false)."); end
    if params.save
        fprintf('Saving outputs set to ON. \n')
    else
        fprintf('Saving outputs set to OFF. \n')
    end
else
    params.save = true; % save output by default
    disp("Saving outputs not defined. Set to ON by default. If you wish to turn it OFF, set input 'save' to false");
end

% Parallel computing
idx = find(strcmpi(varargin,'parpool'));
ps = parallel.Settings;
if ~isempty(idx)
    params.parpool = varargin{idx+1};
    if ~islogical(params.parpool), error("The 'parpool' input should be logical (true or false)."); end
    if params.parpool
        fprintf('Parallel computing set to ON. \n')
        params.parpool = true;
        ps.Pool.AutoCreate = true;
    else
        fprintf('Parallel computing set to OFF. \n')
        params.parpool = false;
        ps.Pool.AutoCreate = false;  % prevents parfor loops from launching parpool mode
    end
else
    params.parpool = false;
    ps.Pool.AutoCreate = false;  % prevents parfor loops from launching parpool mode
    fprintf('Parallel computing not defined: set to OFF by default. \n')
end

% GPU computing 
idx = find(strcmpi(varargin,'gpu'));
if ~isempty(idx)
    params.gpu = varargin{idx+1};
    if ~islogical(params.gpu), error("The 'gpu' input should be logical (true or false)."); end
    fprintf('GPU computing: set to ON. \n')
else
    params.gpu = false;
    fprintf('GPU computing not defined: set to OFF by default. \n')
end
