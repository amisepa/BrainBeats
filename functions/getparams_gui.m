% GETPARAMS_GUI - Get BrainBeats parameters from the user through two windows.
%
% Usage:
%   [params, abort, EEG] = getparams_gui(EEG)
%
% The main window (BRAINBEATS_MAIN_WINDOW, with the BrainBeats logo) shows
% the dataset (or loads one) and selects the analysis ('hep', 'features' or
% 'coherence'), the heart signal and channel(s), for HEP whether
% the measures are computed on the scalp channels, the components or both, and the
% plots and saving. The parameters window that follows sets, depending on
% the analysis: the heartbeat detection (ECG: R-peak threshold, refractory
% period, searchback; PPG: pulse fiducial and, for HEP, the pulse arrival
% time), the EEG preprocessing (artifact removal with ASR/bad epochs + ICA
% or GEDAI; reference: average, infinity, surface Laplacian or none;
% filters; bad channels; ICA method), the HEP options (epoch window,
% baseline, ROI, HRSP/HRPC of all channels, surrogate control, heart channel
% kept), and the features. GUI values are converted to
% the formats expected by BRAINBEATS_PROCESS (see its help for each option).
% Internal, called by brainbeats_process.
%
% Inputs:
%   EEG    - EEGLAB EEG structure (can be empty: loaded in the main window)
% Outputs:
%   params - BrainBeats parameters
%   abort  - true if the user closed or cancelled a window
%   EEG    - the dataset (loaded in the main window if the input was empty)
%
% Copyright (C) - Cedric Cannard, 2023

function [params, abort, EEG] = getparams_gui(EEG)

%% Main window
[params, EEG, abort] = brainbeats_main_window(EEG);
if abort, return; end
abort = true;
isHEP = strcmp(params.analysis,'hep');
isECG = strcmp(params.heart_signal,'ecg');

%% Parameters window
% Popup lists (the first item is the default unless 'value' says otherwise)
cleanMethods = {'ASR / bad epochs + ICA (ICLabel)' 'GEDAI (plugin installed if needed)'};
refMethods = {'Average' 'Infinity (REST)' 'Surface Laplacian (CSD)' 'None'};
lineFreqs = {'50 Hz' '60 Hz'};
filtTypes = {'Zero-phase (default)' 'Causal (minimum-phase)'};
icaMethods = {'Picard (fast)' 'Infomax (default)' 'Infomax, replicable (slow)'};
badEpochs = {'Aggressive (quartiles)' 'Moderate (Grubbs, default)' 'Lax (mean)'};
ppgModes = {'Pulse onsets (valleys, default)' 'Systolic peaks'};
patModes = {'Automatic (recommended)' 'Fixed delay (ms):' 'None'};
hepWins = {'Fixed (default)' 'Adaptive to the heart rate (within-subject only)'};
hepBls = {'None (recommended)' 'Regression (Alday, 2019)'};
spectypes = {'Normalized Lomb-Scargle periodogram (default)' 'Lomb-Scargle periodogram' 'Welch (requires resampling)' 'FFT (requires resampling)'};
hrvnorm = {'Yes' 'No (default)'};
freqrange = '[1 40]';
wintype = {'hamming' 'hann' 'rectwin' 'blackmanharris'};
freqbounds = {'Conventional (e.g., alpha = 8-13 Hz)' 'Individualized (e.g., alpha = 7.8-12.3 Hz)'};
winlen = '2';
eegnorm = {'None (uV^2/Hz)' 'Decibels (default)' 'Decibels + divided by total power'};
asynorm = {'None' 'Divided by total power'};
roiDefault = 'F1 F2 F3 F4 Fz FC1 FC2 FC3 FC4 FCz C1 C2 C3 C4 Cz';
heartRemovals = {'ICA (ICLabel heart components)' 'ECG regression (lags -20 to 20 ms)' 'None'};
heartRemovalKeys = {'ica' 'ecg_regression' 'none'};
if ~isECG   % the regression needs the ECG
    heartRemovals(2) = []; heartRemovalKeys(2) = [];
end

% Callbacks: enable the options of a checked box (userdata = its tag)
cbEEG = "if get(gcbo,'value'), set(findobj(gcbf,'userdata','clean_eeg'),'enable','on'); else, set(findobj(gcbf,'userdata','clean_eeg'),'enable','off'); end";
cbHRV = "if get(gcbo,'value'), set(findobj(gcbf,'userdata','hrv_features'),'enable','on'); else, set(findobj(gcbf,'userdata','hrv_features'),'enable','off'); end";
cbEEGf = "if get(gcbo,'value'), set(findobj(gcbf,'userdata','eeg_features'),'enable','on'); else, set(findobj(gcbf,'userdata','eeg_features'),'enable','off'); end";
cbROI = "tmpEEG = get(gcbf,'userdata'); [~,tmpv] = pop_chansel({tmpEEG.chanlocs.labels},'withindex','on'); if ~isempty(tmpv), set(findobj(gcbf,'tag','hep_roi'),'string',tmpv); end; clear tmpEEG tmpv";
L = {}; G = {};   % uilist and geometry, filled block by block
row3 = [.02 .45 .35];
row4 = [.02 .45 .2 .15];

% Heartbeat detection
if isECG
    L = [L {{'style' 'text' 'string' 'Heartbeat detection (ECG)' 'fontweight' 'bold'}} ...
        {{} {'style' 'text' 'string' 'R-peak threshold (lower = more sensitive):'} {'style' 'edit' 'string' '0.35' 'tag' 'ecg_peakthresh'}} ...
        {{} {'style' 'text' 'string' 'Refractory period (s):'} {'style' 'edit' 'string' '0.25' 'tag' 'ecg_refperiod'}} ...
        {{} {'style' 'text' 'string' 'Search again long gaps at half threshold:'} {'style' 'popupmenu' 'string' {'Yes (default)' 'No'} 'tag' 'ecg_searchback'}}];
    G = [G {1 row3 row3 row3}];
else
    L = [L {{'style' 'text' 'string' 'Heartbeat detection (PPG)' 'fontweight' 'bold'}} ...
        {{} {'style' 'text' 'string' 'Pulse fiducial:'} {'style' 'popupmenu' 'string' ppgModes 'tag' 'ppg_detect_mode'}}];
    G = [G {1 row3}];
    if isHEP
        labs = {EEG.chanlocs.labels};
        ecgLab = labs(contains(lower(labs), {'ecg' 'ekg'}));
        if ~isempty(ecgLab)
            patInfo = sprintf('Automatic: measured from ECG channel %s', ecgLab{1});
        else
            patInfo = 'Automatic: from the EEG cardiac field artifact (no ECG), else 250 ms';
        end
        L = [L {{} {'style' 'text' 'string' 'Pulse arrival time (shifts the pulses to the heartbeats):'} ...
            {'style' 'popupmenu' 'string' patModes 'tag' 'ppg_transit_mode'} {'style' 'edit' 'string' '250' 'tag' 'ppg_transit'}} ...
            {{} {'style' 'text' 'string' patInfo 'fontangle' 'italic'}}];
        G = [G {row4 [.02 .8]}];
    end
end
L = [L {{}}]; G = [G {1}];

% EEG preprocessing
hp = '1'; lp = '40';
if isHEP, hp = '0.5'; lp = '30'; end
L = [L {{'style' 'checkbox' 'string' 'Preprocess the EEG' 'fontweight' 'bold' 'tag' 'clean_eeg' 'callback' cbEEG 'value' 1}}];
G = [G {1}];
L = [L {{} {'style' 'text' 'string' 'Artifact removal:'} {'style' 'popupmenu' 'string' cleanMethods 'tag' 'clean_method' 'userdata' 'clean_eeg'}}];
G = [G {row3}];
L = [L {{} {'style' 'text' 'string' 'Reference:'} {'style' 'popupmenu' 'string' refMethods 'tag' 'ref' 'userdata' 'clean_eeg'}} ...
    {{} {'style' 'text' 'string' 'Power line frequency:'} {'style' 'popupmenu' 'string' lineFreqs 'tag' 'linenoise' 'userdata' 'clean_eeg'}} ...
    {{} {'style' 'text' 'string' 'High-pass filter (Hz; ICA uses a 1-Hz copy):'} {'style' 'edit' 'string' hp 'tag' 'highpass' 'userdata' 'clean_eeg'}} ...
    {{} {'style' 'text' 'string' 'Low-pass filter (Hz):'} {'style' 'edit' 'string' lp 'tag' 'lowpass' 'userdata' 'clean_eeg'}} ...
    {{} {'style' 'text' 'string' 'Filter type:'} {'style' 'popupmenu' 'string' filtTypes 'tag' 'filttype' 'userdata' 'clean_eeg'}} ...
    {{} {'style' 'text' 'string' 'Bad channels: minimum correlation (.5-.9):'} {'style' 'edit' 'string' '.65' 'tag' 'corrThresh' 'userdata' 'clean_eeg'}}];
G = [G {row3 row3 row3 row3 row3 row3}];
if isHEP
    L = [L {{} {'style' 'text' 'string' 'Bad epochs detection:'} {'style' 'popupmenu' 'string' badEpochs 'tag' 'detectMethod' 'value' 2 'userdata' 'clean_eeg'}} ...
        {{} {'style' 'text' 'string' 'Heart artifact removal (ICA: min. ICLabel heart probability, %):'} ...
            {'style' 'popupmenu' 'string' heartRemovals 'tag' 'heart_removal' 'userdata' 'clean_eeg'} {'style' 'edit' 'string' '75' 'tag' 'conf_thresh' 'userdata' 'clean_eeg'}}];
    G = [G {row3 row4}];
else
    L = [L {{} {'style' 'text' 'string' 'ASR threshold for bad segments (SD):'} {'style' 'edit' 'string' '30' 'tag' 'asr_cutoff' 'userdata' 'clean_eeg'}}];
    G = [G {row3}];
end
L = [L {{} {'style' 'text' 'string' 'ICA method:'} {'style' 'popupmenu' 'string' icaMethods 'tag' 'icamethod' 'value' 2 'userdata' 'clean_eeg'}} {{}}];
G = [G {row3 1}];

% HEP
if isHEP
    L = [L {{'style' 'text' 'string' 'HEP' 'fontweight' 'bold'}} ...
        {{} {'style' 'text' 'string' 'Epoch window (ms):'} {'style' 'popupmenu' 'string' hepWins 'tag' 'hep_window_mode'} {'style' 'edit' 'string' '[-300 600]' 'tag' 'hep_window'}} ...
        {{} {'style' 'text' 'string' 'Baseline correction (window in ms):'} {'style' 'popupmenu' 'string' hepBls 'tag' 'hep_baseline'} {'style' 'edit' 'string' '[-150 -50]' 'tag' 'hep_baseline_win'}}];
    G = [G {1 row4 row4}];
    if ~isfield(params,'hep_level') || ~strcmp(params.hep_level,'ics')
        L = [L {{} {'style' 'text' 'string' 'ROI channels (plots):'} {'style' 'edit' 'string' roiDefault 'tag' 'hep_roi'} {'style' 'pushbutton' 'string' '...' 'callback' cbROI}}];
        G = [G {row4}];
    end
    L = [L {{} {'style' 'text' 'string' 'HRSP and HRPC frequency range (Hz):'} {'style' 'edit' 'string' '[4 30]' 'tag' 'hep_tf_freqs'}} ...
        {{} {'style' 'text' 'string' 'Surrogate heartbeat control (number of surrogates, 0 = none):'} {'style' 'edit' 'string' '100' 'tag' 'hep_surrogates'}} ...
        {{} {'style' 'checkbox' 'string' 'Keep the heart channel in the output (its own units)' 'tag' 'keep_heart' 'value' 0}}];
    G = [G {row3 row3 [.02 .8]}];
end

% Features
if strcmp(params.analysis,'features')
    L = [L {{'style' 'checkbox' 'string' 'HRV features' 'fontweight' 'bold' 'tag' 'hrv_features' 'callback' cbHRV 'value' 1}} ...
        {{} {'style' 'checkbox' 'string' 'Time domain (SDNN, RMSSD, pNN50)' 'tag' 'hrv_time' 'value' 1 'userdata' 'hrv_features'}} ...
        {{} {'style' 'checkbox' 'string' 'Frequency domain (ULF, VLF, LF, HF, LF/HF, total power)' 'tag' 'hrv_frequency' 'value' 1 'userdata' 'hrv_features'} ...
            {'style' 'edit' 'string' '' 'tag' 'hrv_freq_opts' 'visible' 'off'} {'style' 'pushbutton' 'string' 'Options' 'callback' {@hrvfreqparam spectypes hrvnorm} 'userdata' 'hrv_features'}} ...
        {{} {'style' 'checkbox' 'string' 'Nonlinear domain (Poincare, PRSA, entropy, fractal dimension)' 'tag' 'hrv_nonlinear' 'value' 1 'userdata' 'hrv_features'}} ...
        {{'style' 'checkbox' 'string' 'EEG features' 'fontweight' 'bold' 'tag' 'eeg_features' 'callback' cbEEGf 'value' 1}} ...
        {{} {'style' 'checkbox' 'string' 'Time domain (RMS, variance, skewness, kurtosis, IQR)' 'tag' 'eeg_time' 'value' 1 'userdata' 'eeg_features'}} ...
        {{} {'style' 'checkbox' 'string' 'Frequency domain (band power, IAF, asymmetry)' 'tag' 'eeg_frequency' 'value' 1 'userdata' 'eeg_features'} ...
            {'style' 'edit' 'string' '' 'tag' 'eeg_freq_opts' 'visible' 'off'} {'style' 'pushbutton' 'string' 'Options' 'callback' {@eegfreqparam freqrange wintype freqbounds winlen eegnorm asynorm} 'userdata' 'eeg_features'}} ...
        {{} {'style' 'checkbox' 'string' 'Nonlinear domain (fuzzy entropy, fractal dimension; slow)' 'tag' 'eeg_nonlinear' 'value' 1 'userdata' 'eeg_features'}} ...
        {{'style' 'checkbox' 'string' 'Use parallel computing' 'tag' 'parpool' 'value' 0}} ...
        {{'style' 'checkbox' 'string' 'Use GPU computing' 'tag' 'gpu' 'value' 0}}];
    G = [G {1 [.02 .8] [.02 .6 .01 .15] [.02 .8] 1 [.02 .8] [.02 .6 .01 .15] [.02 .8] 1 1}];
end

uilist = L;     % one cell per control ({} = empty slot)
uigeom = G;     % one vector per row
titles = struct('hep','HEP, HRSP and HRPC','features','features','coherence','brain-heart coherence');
[res,~,~,p2] = inputgui(uigeom, uilist, 'pophelp(''brainbeats_process'')', ...
    sprintf('BrainBeats: parameters for %s', titles.(params.analysis)), EEG);
if isempty(res), return; end
abort = false;

%% Convert the GUI values to BRAINBEATS_PROCESS parameters
f = fieldnames(p2);
for i = 1:numel(f), params.(f{i}) = p2.(f{i}); end

num = @(v) str2num(v); %#ok<ST2NM>
logicals = {'clean_eeg' 'keep_heart' 'hrv_features' 'hrv_time' 'hrv_frequency' 'hrv_nonlinear' ...
    'eeg_features' 'eeg_time' 'eeg_frequency' 'eeg_nonlinear' 'parpool' 'gpu'};
for i = 1:numel(logicals)
    if isfield(params, logicals{i}), params.(logicals{i}) = logical(params.(logicals{i})); end
end
numbers = {'ecg_peakthresh' 'ecg_refperiod' 'highpass' 'lowpass' 'corrThresh' 'asr_cutoff' 'hep_surrogates'};
for i = 1:numel(numbers)
    if isfield(params, numbers{i}), params.(numbers{i}) = str2double(params.(numbers{i})); end
end
if isfield(params,'ecg_searchback'), params.ecg_searchback = params.ecg_searchback == 1; end
if isfield(params,'ppg_detect_mode')
    modes = {'valleys' 'peaks'}; params.ppg_detect_mode = modes{params.ppg_detect_mode};
end
if isfield(params,'ppg_transit_mode')
    switch params.ppg_transit_mode
        case 1, params.ppg_transit = 'auto';
        case 2, params.ppg_transit = str2double(params.ppg_transit);
        case 3, params.ppg_transit = 'off';
    end
    params = rmfield(params,'ppg_transit_mode');
end
if isfield(params,'clean_method')
    methods = {'asr_ica' 'gedai'}; params.clean_method = methods{params.clean_method};
end
if isfield(params,'heart_removal'), params.heart_removal = heartRemovalKeys{params.heart_removal}; end
if isfield(params,'ref')
    refs = {'average' 'infinity' 'csd' 'off'}; params.ref = refs{params.ref};
end
if isfield(params,'linenoise'), params.linenoise = 10*(4 + params.linenoise); end   % 50 or 60
if isfield(params,'filttype')
    types = {'noncausal' 'causal'}; params.filttype = types{params.filttype};
end
if isfield(params,'detectMethod')
    dm = {'quartiles' 'grubbs' 'mean'}; params.detectMethod = dm{params.detectMethod};
end
if isfield(params,'lowpass') && isfield(params,'linenoise') && params.linenoise < params.lowpass
    warning('The low-pass filter (%g Hz) is above the power line frequency: a notch filter is applied at %g Hz.', params.lowpass, params.linenoise)
end
if isfield(params,'conf_thresh'), params.conf_thresh = str2double(params.conf_thresh)/100; end

% HEP
if isfield(params,'hep_window_mode')
    if params.hep_window_mode == 2
        params.hep_window = 'adaptive';
    else
        params.hep_window = num(params.hep_window);
        if numel(params.hep_window) ~= 2 || params.hep_window(1) >= 0
            error('The HEP epoch window must be [start end] in ms, with start < 0 (e.g. [-300 600]).')
        end
    end
    params = rmfield(params,'hep_window_mode');
end
if isfield(params,'hep_baseline')
    if params.hep_baseline == 2
        params.hep_baseline = 'regression';
        params.hep_baseline_win = num(params.hep_baseline_win);
    else
        params.hep_baseline = 'none';
        params = rmfield(params,'hep_baseline_win');
    end
end
if isfield(params,'hep_roi')
    params.hep_roi = strsplit(strtrim(params.hep_roi));
end
if isfield(params,'hep_tf_freqs'), params.hep_tf_freqs = num(params.hep_tf_freqs); end
if isfield(params,'hep_surrogates') && (isnan(params.hep_surrogates) || params.hep_surrogates < 0)
    error('The number of surrogates must be 0 (none) or more.')
end

% HRV frequency options (from the 'Options' button)
if isfield(params,'hrv_freq_opts')
    o = params.hrv_freq_opts;
    if iscell(o) && ~isempty(o)
        specs = {'LombScargle_norm' 'LombScargle' 'welch' 'fft'};
        params.hrv_spec = specs{str2double(o{find(strcmp(o,'hrvspec'))+1})};
        params.hrv_overlap = str2double(o{find(strcmp(o,'winoverlap'))+1})/100;
        params.hrv_norm = strcmp(o{find(strcmp(o,'hrvnorm'))+1}, '1');
    end
    params = rmfield(params,'hrv_freq_opts');
end

% EEG frequency options (from the 'Options' button)
if isfield(params,'eeg_freq_opts')
    o = params.eeg_freq_opts;
    if iscell(o) && ~isempty(o)
        params.eeg_frange = num(o{find(strcmp(o,'frange'))+1});
        params.eeg_wintype = o{find(strcmp(o,'wintype'))+1};
        params.eeg_winoverlap = str2double(o{find(strcmp(o,'winoverlap'))+1});
        params.eeg_winlen = str2double(o{find(strcmp(o,'winlen'))+1});
        if contains(o{find(strcmp(o,'freqbounds'))+1}, 'Conventional')
            params.eeg_freqbounds = 'conventional';
        else
            params.eeg_freqbounds = 'individualized';
        end
        tmp = o{find(strcmp(o,'eegnorm'))+1};
        if contains(lower(tmp),'none'), params.eeg_norm = 0;
        elseif contains(tmp,'default'), params.eeg_norm = 1;
        else, params.eeg_norm = 2;
        end
        params.asy_norm = ~contains(o{find(strcmp(o,'asynorm'))+1}, 'None');
    end
    params = rmfield(params,'eeg_freq_opts');
end

% Domains of a disabled feature group are off (as in getparams_cmd)
if isfield(params,'hrv_features') && ~params.hrv_features
    [params.hrv_time, params.hrv_frequency, params.hrv_nonlinear] = deal(false);
end
if isfield(params,'eeg_features') && ~params.eeg_features
    [params.eeg_time, params.eeg_frequency, params.eeg_nonlinear] = deal(false);
end

% EEG preprocessing off: drop its options
if ~params.clean_eeg
    drop = {'clean_method' 'ref' 'linenoise' 'highpass' 'lowpass' 'filttype' 'corrThresh' 'detectMethod' 'asr_cutoff' 'icamethod' 'heart_removal' 'conf_thresh'};
    params = rmfield(params, drop(isfield(params, drop)));
end
if ~isfield(params,'parpool'), params.parpool = false; end
if ~isfield(params,'gpu'), params.gpu = false; end
params.heart = true;
params.eeg = true;


%% Button to get options for HRV frequency estimation

function hrvfreqparam(~, ~, spectypes, hrvnorm)
uilist = {
    {'style' 'text' 'string' 'Method to estimate HRV power' } {'style' 'popupmenu' 'string' spectypes 'tag' 'hrv_spec'}  ...
    {} ...
    {'style' 'text' 'string' 'Window overlap (in %):' } {'style' 'edit' 'string' '25' 'tag' 'hrv_overlap' }  ...
    {} ...
    {'style' 'text' 'string' 'Normalize power' } {'style' 'popupmenu' 'string' hrvnorm 'tag' 'hrv_normc' 'value' 2}  ...
    };
uigeom = { [.2 .3] .3 [.2 .1] .3 [.2 .2] };
result = inputgui(uigeom, uilist, 'help(''get_hrv_features'')', 'HRV frequency domain parameters');
if isempty(result), return, end
% Store the choices in the hidden 'hrv_freq_opts' field of the parent window
out = { 'hrvspec' num2str(result{1}) 'winoverlap' result{2} 'hrvnorm' num2str(result{3}) };
set(findobj(gcbf, 'tag', 'hrv_freq_opts'), 'string', out );


%% Button to get options for EEG frequency features

function eegfreqparam(~, ~, freqrange, wintype, freqbounds, winlen, eegnorm, asynorm)
uilist = {
    {'style' 'text' 'string' 'Overall frequency range' 'fontweight' 'bold'} {'style' 'edit' 'string' freqrange 'tag' 'eeg_freqrange'}  ...
    {} ...
    {'style' 'text' 'string' 'Window type:' 'fontweight' 'bold'} {'style' 'popupmenu' 'string' wintype 'tag' 'eeg_wintype' }  ...
    {} ...
    {'style' 'text' 'string' 'Window overlap (in %):' 'fontweight' 'bold'} {'style' 'edit' 'string' '50' 'tag' 'eeg_overlap' }  ...
    {} ...
    {'style' 'text' 'string' 'Window length (in s):' 'fontweight' 'bold'} {'style' 'edit' 'string' winlen 'tag' 'eeg_winlen' }  ...
    {} ...
    {'style' 'text' 'string' 'Frequency bands' 'fontweight' 'bold'} {'style' 'popupmenu' 'string' freqbounds 'tag' 'eeg_freqbounds' }  ...
    {} ...
    {'style' 'text' 'string' 'Band-power normalization' 'fontweight' 'bold'} {'style' 'popupmenu' 'string' eegnorm 'tag' 'eeg_norm' 'value' 2}  ...
    {} ...
    {'style' 'text' 'string' 'Alpha asymmetry normalization' 'fontweight' 'bold'} {'style' 'popupmenu' 'string' asynorm 'tag' 'asy_norm' }  ...
    };
uigeom = { [.1 .1] .3 [.1 .1] .3 [.1 .1] .3 [.1 .1] .3 [.2 .2] .3 [.2 .2] .3 [.1 .1] };
result = inputgui(uigeom, uilist, 'help(''get_eeg_features'')', 'EEG frequency domain parameters');
if isempty(result), return, end
out = { 'frange' result{1} 'wintype' wintype{result{2}} 'winoverlap' result{3} ...
    'winlen' result{4} 'freqbounds' freqbounds{result{5}} 'eegnorm' eegnorm{result{6}} ...
    'asynorm' asynorm{result{7}} };
set(findobj(gcbf, 'tag', 'eeg_freq_opts'), 'string', out );  % parsed in the main function
