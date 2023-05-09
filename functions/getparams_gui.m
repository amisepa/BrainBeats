% Get parameters specified by user via GUI.
% 
% Cedric Cannard, April 2023 

function params = getparams_gui(EEG)

% dropdown options to select
dataTypes = { 'ECG' 'PPG' };
rrArtMethod = { 'Shape-preserving piecewise cubic interpolation (default)' ...
    'Linear interpolation' 'Cubic interpolation' 'Nearest neighbor interpolation' ...
    'Next neighbor interpolation' 'Previous interpolation' 'Spline interpolation'  ...
    'Cubic convolution interpolation' 'Modified Akima cubic interpolation' 'Remove'};
analysisTypes = { 'Heartbeat-evoked potentials (HEP)' 'Extract EEG & HRV features from continuous data' 'Remove heart artifacts from EEG signals'};
cleanEEG = {'Yes' 'No (already processed)'};

% callback functions to allow features to appear
slctAnalysis = "if get(gcbo,'value') <= 2, set(findobj(gcbf,'userdata','analysis'),'enable','on'); else, set(findobj(gcbf,'userdata','analysis'),'enable','off'); end";
slctChan = "tmpEEG = get(gcbf, 'userdata'); tmpchanlocs = tmpEEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels},'withindex','on'); set(findobj(gcbf,'tag','heart_channels'),'string',tmpval); clear tmp tmpEEG tmpchanlocs tmpval";
computeEEG = "if get(gcbo,'value'), set(findobj(gcbf,'userdata','eeg'),'enable','on'); else, set(findobj(gcbf,'userdata','eeg'),'enable','off'); end";
computeHRV = "if get(gcbo,'value'), set(findobj(gcbf,'userdata','hrv'),'enable','on'); else, set(findobj(gcbf,'userdata','hrv'),'enable','off'); end";
uigeom = { [.3 .6] 1 [.6 .3] 1 [.3 .3 .3] 1 [.3 .6] 1 [.6 .3] 1 ...
    [.15 .6] [.2 .6] [.2 .6] 1 ...
    [.15 .6] [.2 .6] [.2 .6] [.2 .6] ...
    [1 1]};

uilist = { ...
    {'style' 'text' 'string' 'Analysis to run:' 'fontweight' 'bold' } {'style' 'popupmenu' 'string' analysisTypes 'tag' 'analysis' 'callback' slctAnalysis} ...
    {} ...
    {'style' 'text' 'string' 'Heart data type' 'fontweight' 'bold'} {'style' 'popupmenu' 'string' dataTypes 'tag' 'heart_signal' 'value' 1} ...
    {} ...
    {'style' 'text' 'string' 'Select ECG/PPG channel(s):' 'fontweight' 'bold'} {'style' 'edit' 'tag' 'heart_channels'} {'style' 'pushbutton' 'string'  '...', 'enable' 'on' 'callback'  slctChan } ...
    {} ...
    {'style' 'text' 'string' 'What to do with RR artifacts:' 'fontweight' 'bold' } {'style' 'popupmenu' 'string' rrArtMethod 'userdata' 'analysis' 'enable' 'off' 'value' 1 'tag' 'rr_correct'} ...
    {} ...
    {'style' 'text' 'string' 'Clean EEG data?' 'fontweight' 'bold' } {'style' 'popupmenu' 'string' cleanEEG 'tag' 'clean_eeg'} ...
    {} ...
    {} {'style' 'checkbox' 'string' ' EEG features' 'fontweight' 'bold' 'tag' 'eeg' 'userdata' 'analysis' 'enable' 'off' 'callback' computeEEG 'value' 0} ...
    ...
    {} {'style' 'checkbox' 'tag' 'eeg_frequency' 'string' ' Frequency domain (Delta, Theta, Alpha, Beta, IAF, Asymmetry, Coherence)' 'value' 1 'userdata' 'eeg' 'enable' 'off'} ...
    ...
    {} {'style' 'checkbox' 'tag' 'eeg_nonlinear' 'string' ' Nonlinear domain (entropy)' 'value' 1 'userdata' 'eeg' 'enable' 'off'} ...
    {} ...
    {} {'style' 'checkbox' 'string' ' HRV features' 'fontweight' 'bold' 'tag' 'hrv' 'userdata' 'analysis' 'enable' 'off' 'callback' computeHRV 'value' 0} ...
    ...
    {} {'style' 'checkbox' 'tag' 'hrv_time' 'string' ' Time domain (NN features, SDNN, RMSSD, pNN50)' 'value' 1 'userdata' 'hrv' 'enable' 'off'} ...
    ...
    {} {'style' 'checkbox' 'tag' 'hrv_frequency' 'string' ' Frequency domain (ULF, VLF, LF, HF, LF/HF, Total power)' 'value' 1 'userdata' 'hrv' 'enable' 'off'} ...
    ...
    {} {'style' 'checkbox' 'tag' 'hrv_nonlinear' 'string' ' Nonlinear domain (Poincaré, PRSA, entropy)' 'value' 1 'userdata' 'hrv' 'enable' 'off'} ...
    ...
    {'style' 'checkbox' 'string' ' Plot outputs' 'fontweight' 'bold' 'tag' 'vis' 'value' 1} {} ...
    };

% Launch GUI and get parameters from user
[res,~,~,params] = inputgui(uigeom,uilist,'pophelp(''pop_BrainBeats'')','BrainBeats EEGLAB plugin',EEG);

% Exit if no input
if sum([res{:}]) == 0, return; end

% Analysis choice and check data compatibility for that analysis
if params.analysis == 1
    params.analysis = 'hep';
elseif params.analysis == 2
    params.analysis = 'features';
elseif params.analysis == 3
    params.analysis = 'rm_heart';
else
    % If user didn't select any analysis, select based on data format
    if length(size(EEG.data)) == 2
        params.analysis = 'features';
        disp("You did not select which analysis you want to perform. Continuous data detected: performing feature-based analysis.")
    elseif length(size(EEG.data)) == 3
        params.analysis = 'hep';
        disp("You did not select which analysis you want to perform. Epoched data detected: performing HEP analysis.")
    end
end

% Extract heart signal type
params.heart_signal = lower(dataTypes{params.heart_signal});

% Extract heart channel names
if ~isempty(params.heart_channels)
    if contains(params.heart_channels,' ')
        params.heart_channels = split(params.heart_channels);
    end
else
    error('You must select your ECG/PPG channels from the list')
end

% RR artfact correction method
if ~isempty(params.rr_correct)
    switch params.rr_correct
        case 1
            params.rr_correct = 'pchip';
        case 2            
            params.rr_correct = 'linear';
        case 3
            params.rr_correct = 'cubic';
        case 4
            params.rr_correct = 'nearest';
        case 5
            params.rr_correct = 'next';
        case 6
            params.rr_correct = 'previous';
        case 7
            params.rr_correct = 'spline';
        case 8
            params.rr_correct = 'cubic';
        case 9
            params.rr_correct = 'makima';
        case 'Remove'
            params.rr_correct = 'remove';
    end
end

% Clean EEG?
if params.clean_eeg == 1
    params.clean_eeg = true;
else
    params.clean_eeg = false;
end
    
% Visualize outputs
params.vis = logical(params.vis);

end
