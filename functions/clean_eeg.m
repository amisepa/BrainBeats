% CLEAN_EEG - Preprocess EEG data in two stages, selected by params.clean_eeg_step.
%
% Usage:
%   [EEG, params] = clean_eeg(EEG, params)
%
% Stage 0 (params.clean_eeg_step == 0): filter, re-reference, remove bad channels
%   - High-pass (params.highpass, default 0.5 Hz for 'hep' to keep the slow HEP
%     components, 1 Hz otherwise) and low-pass (params.lowpass, default 30 Hz)
%     FIR filters (pop_eegfiltnew). Zero-phase by default; minimum-phase
%     (causal) if params.filttype = 'causal'.
%   - Re-reference (params.ref: 'average' (default, full-rank CAR via apply_car),
%     'infinity' (REST), 'csd' or 'off'). Skipped with < 30 channels. With
%     'csd', the data are average-referenced here and the surface Laplacian
%     is applied after artifact removal by the caller (ICLabel needs
%     average-referenced data; see CSD_TRANSFORM). If 'infinity' fails, CAR
%     is used instead.
%   - Remove flat channels (params.flatline, default 5 s) and bad channels with
%     clean_channels (params.corrThresh = .65, line noise 15 SD, params.maxBad
%     = .33, 100 RANSAC samples). The location-free clean_channels_nolocs is
%     used for MEG labels or when channel locations are unusable.
%   - Band-stop at params.linenoise +/- 3 Hz if it is below the low-pass.
%   - params.clean_method = 'gedai': GEDAI (Ros et al., 2025; plugin installed
%     by RUN_CHECKS) removes artifacts from the continuous data here, instead
%     of ASR/bad epochs + ICA in stage 1 (bad HEP epochs are still removed).
%     'infinity' is then applied after GEDAI (which re-references to average).
% Stage 1 (params.clean_eeg_step == 1): remove artifacts, interpolate, ICA
%   - 'hep': remove bad epochs with find_badTrials (params.detectMethod,
%     default 'grubbs'). 'features', 'rm_heart', 'coherence': remove bad
%     segments with ASR (params.asr_cutoff, default 50 SD; params.asr_mem,
%     default .85 of available RAM).
%   - Interpolate removed channels back to params.orichanlocs (> 10 channels).
%   - ICA at the effective data rank to avoid ghost ICs (Kim et al. 2023).
%     With a high-pass below 1 Hz, ICA is fitted on a 1-Hz high-passed copy
%     (params.ica_source for epoched data, prepared by RUN_HEP) and the
%     weights are applied to the data: slow drifts degrade ICA and ICLabel
%     (Winkler et al., 2015). With 'gedai', ICA runs only for the IC-level
%     HEP measures (params.hep_level 'ics' or 'both'), and no component is removed.
%     params.icamethod: 1 = Picard, 2 = extended Infomax (default),
%     3 = extended Infomax with lrate 1e-5 and maxsteps 2000 (slower, more
%     replicable).
%   - ICLabel, then remove muscle, eye, heart, line and channel-noise ICs:
%       'rm_heart': .99/.90/-/.99/.99 (heart left for remove_heartcomp)
%       'hep':      .99/.90/.75/.99/.99 (heart: params.conf_thresh; none if
%                   params.heart_removal is 'ecg_regression' or 'none')
%       otherwise:  .95/.95/.99/.99/.99
%
% Inputs:
%   EEG    - EEGLAB EEG structure (EEG channels only)
%   params - BrainBeats parameters (also reads analysis, vis_cleaning, gpu)
% Outputs:
%   EEG    - cleaned EEG structure
%   params - with defaults filled in, clean_eeg_step incremented, and what was
%            removed: bad_channels (logical mask), removed_eeg_channels (labels),
%            removed_eeg_trials, removed_eeg_segments ([start end] samples),
%            removed_eeg_components
%
% Copyright (C) - Cedric Cannard, 2023

function [EEG, params] = clean_eeg(EEG, params)

% General parameters
if isfield(params,'ref')
    reref = params.ref;
else
    reref = 'average'; % 'average' (default), 'infinity', 'csd', 'off'
    params.ref = 'average';    % for user
end
if isfield(params,'highpass') && ~isempty(params.highpass)
    highpass = params.highpass;
elseif strcmp(params.analysis,'hep')
    highpass = 0.5; % HEP: keeps the slow HEP components (ICA uses a 1-Hz copy)
    params.highpass = 0.5;
else
    highpass = 1; % default = 1 Hz
    params.highpass = 1;    % for user
end
if isfield(params,'clean_method') && strcmpi(params.clean_method,'gedai')
    useGEDAI = true;
else
    useGEDAI = false;
    params.clean_method = 'asr_ica';
end
if isfield(params,'lowpass')
    lowpass = params.lowpass;
else
    lowpass = 30;   % default = 30 Hz: removes 50/60 Hz line noise with a low
                    % filter order (fewer ripples, faster), speeds up ICA
                    % and smooths ERPs
    params.lowpass = 30;    % for user
end
if isfield(params,'filttype')
    if strcmpi(params.filttype,'causal')
        causalfilt = true;  % causal (nonlinear) minimum-phase filter
    else 
        causalfilt = false; % zero-phase (linear) noncausal filter
    end
else
    % default filter
    causalfilt = false; % zero-phase noncausal filter
    params.causalfilt = 0;    % to export for user knowledge
end
if isfield(params,'gpu')
    usegpu = params.gpu;
else
    usegpu = false; % default
    params.gpu = 0;    % for user
end

% Channel removal parameters
if isfield(params,'flatline')
    flatline = params.flatline;
else
    flatline = 5; % max flat segment to remove channel (default = 5 s)
    params.flatline = 5;    % for user
end
if isfield(params,'corrThresh')
    corrThresh = params.corrThresh;
else
    corrThresh = .65;   % correlation threshold to be considered bad (default = .65)
    params.corrThresh = .65;    % for user
end
if isfield(params,'maxBad')
    maxBad = params.maxBad;
else
    maxBad = .33;       % max tolerated portion of channel to be bad before removal (default = .33)
    params.maxBad = .33;    % for user
end

% HEP parameters to remove bad epochs
if isfield(params,'detectMethod')
    detectMethod = params.detectMethod;
else
    detectMethod = 'grubbs';   % isoutlier method: 'median' (more aggressive), 'grubbs' (moderate; default), 'mean' (more lax)
end

% ASR parameters
if isfield(params,'asr_cutoff')
    asr_cutoff = params.asr_cutoff;
else
    asr_cutoff = 50;  % main ASR SD cutoff (lower = more aggressive, higher = more lax)
end
if isfield(params,'asr_mem')
    asr_mem = params.asr_mem;
else
    asr_mem = .85;     % available RAM to use for ASR (.85 = 85% of available RAM)
end

% ICA parameters
if isfield(params,'icamethod')
    icamethod = params.icamethod;
else
    icamethod = 2;  % 1 = Picard (fast), 2 = extended Infomax, 3 = replicable Infomax (slowest)
end

% Filter, re-reference, and remove bad channels
if params.clean_eeg_step == 0
    
    % High-pass filter (removes slow drifts)
    EEG = pop_eegfiltnew(EEG,'locutoff',highpass,'minphase',causalfilt);

    % Low-pass filter
    EEG = pop_eegfiltnew(EEG,'hicutoff',lowpass,'minphase',causalfilt);

    % Re-reference (average, infinity/REST or CSD)
    % Candia-Rivera, Catrambone, & Valenza (2021). The role of EEG reference 
    % in the assessment of functional brain–heart interplay: From 
    % methodology to user guidelines. Journal of Neuroscience Methods.
    if ~strcmp(reref,'off')
        if EEG.nbchan < 30
            warndlg('Cannot reference these EEG data to infinity or average or Surface Laplacian, not validated with less than 30 channels.')
            warning('Cannot reference these EEG data to infinity or average or Surface Laplacian, not validated with less than 30 channels.')
        else
            if strcmp(reref,'infinity') && useGEDAI
                fprintf('Re-referencing to infinity after GEDAI. \n')
            elseif strcmp(reref,'infinity')
                fprintf('Re-referencing EEG data to infinity. \n')
                try
                    EEG = rest_ref(EEG);
                    params.ref_applied = 'rest';
                catch ME
                    warning('Re-reference to infinity failed (%s): using the common average reference instead. Please report it: https://github.com/amisepa/BrainBeats/issues', ME.message)
                    EEG = apply_car(EEG);  % preserving effective data rank
                end
            elseif strcmp(reref,'average')
                fprintf('Re-referencing EEG data to average. \n')
                EEG = apply_car(EEG);  % preserving effective data rank
            elseif strcmp(reref,'csd')
                fprintf('Re-referencing EEG data to average (the surface Laplacian is applied after artifact removal). \n')
                EEG = apply_car(EEG);
            end
        end
    end
    
    % Remove bad channels
    win_length = [];     % window length to scan channels ([] = clean_channels default)
    line_thresh = 15;    % line noise threshold in SD (default = 15)
    nSamp = 100;        % number of RANSAC samples (~50-500; higher is slower but more accurate and replicable)
    EEG.etc.clean_channel_mask = true(1,EEG.nbchan);
    oriEEG = EEG;
    EEG = clean_flatlines(EEG,flatline);   % remove channels that have flat lines
    try 
        if any(contains(lower({EEG.chanlocs.labels}), 'meg'))
            warning("MEG data detected. Defaulting to 2nd bad channel detection method that does not leverage EEG electrode locations.")
            EEG = clean_channels_nolocs(EEG,0.45,0.1,win_length,.4);
        else
            EEG = clean_channels(EEG,corrThresh,line_thresh,win_length,maxBad,nSamp); 
        end
    catch
        warning('Your dataset has incorrect electrode locations. Using the location-free algorithm to remove bad EEG channels.');
        EEG = clean_channels_nolocs(EEG,0.45,0.1,win_length,.4);
    end
    badChan = ~ismissing({oriEEG.chanlocs.labels}, {EEG.chanlocs.labels});
    params.bad_channels = badChan; % will be used for feature outputs
    fprintf(1, 'Bad EEG channels removed from data: ');
    fprintf(1, '%s ', oriEEG.chanlocs(badChan).labels )
    fprintf(1, '\n')
    EEG.etc.clean_channel_mask(badChan) = false; 
    badChan = { oriEEG.chanlocs(badChan).labels };
    EEG = pop_select(EEG,'nochannel', badChan);
    
    % Store in params if users want that information
    params.removed_eeg_channels = badChan;

    % Visualize removed channels
    if ~isempty(badChan) && params.vis_cleaning
        try
            vis_artifacts(EEG,oriEEG,'ShowSetname',false);
        catch
            warning('failed to plot bad channels or artifacts with vis_artifacts(). Setting show_events to off and trying againg')
            try
                vis_artifacts(EEG,oriEEG,'ShowSetname',false,'ShowEvents',false);
            catch
                warning("vis_artifacts failed to plot the removed channels or artifacts. Please submit an issue on EEGLAB's page: https://github.com/sccn/eeglab/issues")
            end
        end
        try icadefs; set(gcf, 'color', BACKCOLOR); catch; end     % eeglab background color
        set(gcf,'Toolbar','none','Menu','none');  % remove toolbar and menu
        finish_figure(gcf)
        set(gcf,'Name','EEG channels removed','NumberTitle', 'Off')  % change figure name
    else
        disp("No bad channels detected.")
    end
    
    % Band-stop (notch) filter only if line noise is below the low-pass cutoff
    if isfield(params,'linenoise') && params.linenoise<lowpass
        EEG = pop_eegfiltnew(EEG, 'locutoff',params.linenoise-3, ...
            'hicutoff',params.linenoise+3,'revfilt',1,'filtorder',500);
    end

    % GEDAI on the continuous data (replaces ASR/ICA of stage 1)
    if useGEDAI
        EEG = run_gedai(EEG, params);
        if strcmp(reref,'infinity') && EEG.nbchan >= 30
            try
                EEG = rest_ref(EEG);
                params.ref_applied = 'rest';
            catch ME
                warning('Re-reference to infinity failed (%s): the data stay average-referenced (GEDAI output).', ME.message)
            end
        end
    end

    % update tracker
    params.clean_eeg_step = 1;

% Remove bad epochs (HEP) or bad segments (continuous data), then ICA
elseif params.clean_eeg_step == 1
    
    disp('----------------------------------------------')
    fprintf('              Cleaning EEG data \n')
    disp('----------------------------------------------')

    % HEP (remove bad epochs)
    if strcmp(params.analysis, 'hep')
        
        % Detect and remove bad epochs (and from the copy used for ICA)
        badTrials = find_badTrials(EEG, detectMethod, params.vis_cleaning);
        EEG = pop_rejepoch(EEG, badTrials, 0);
        if isfield(params,'ica_source') && ~isempty(params.ica_source) && ~isempty(badTrials)
            params.ica_source = pop_rejepoch(params.ica_source, badTrials, 0);
        end

        % Store in params if users want that information
        params.removed_eeg_trials = badTrials;
        
    % GEDAI already cleaned the continuous data: no segment removed
    elseif useGEDAI
        params.removed_eeg_segments = [];

    % ASR on continuous data
    elseif contains(params.analysis, {'features' 'rm_heart' 'coherence'})
        
        % Identify artifacts using ASR
        oriEEG = EEG;
        try
            m = memory; maxmem = round(asr_mem*(m.MemAvailableAllArrays/1000000),1);  % asr_mem fraction of available memory (in MB); memory() is Windows-only
            cleanEEG = clean_asr(EEG,asr_cutoff,[],[],[],[],[],[],usegpu,false,maxmem);
        catch
            warning("Failed to use high RAM to run ASR faster. Defaulting back to default values (ASR will just be slower).")
            cleanEEG = clean_asr(EEG,asr_cutoff,[],[],[],[],[],[],usegpu,false,[]);
        end
        
        % Samples modified by ASR = artifacts; convert to [start end] intervals
        mask = sum(abs(EEG.data-cleanEEG.data),1) > 1e-10;
        EEG.etc.clean_sample_mask = true(1, length(mask)); % initialize all samples as clean
        badData = reshape(find(diff([false mask false])), 2, [])';
        badData(:, 2) = badData(:, 2) - 1;
        % keep very short artifacts (< 10 samples) in the data
        if ~isempty(badData)  
            smallIntervals = diff(badData')' < 10;
            badData(smallIntervals, :) = [];
            for i = 1:size(badData, 1)
                EEG.etc.clean_sample_mask(badData(i, 1):badData(i, 2)) = false;
            end
        end

        % Remove them from data
        EEG = pop_select(EEG,'nopoint',badData);
        fprintf('%g %% of data were considered to be artifacts and were removed. \n', (1-EEG.xmax/oriEEG.xmax)*100)
        
        % Store in params if users want that information
        params.removed_eeg_segments = badData;

        % Plot what has been removed
        if params.vis_cleaning
            vis_artifacts(EEG,oriEEG,'ShowSetname',false);
            try icadefs; set(gcf, 'color', BACKCOLOR); catch; end     % eeglab background color
            set(gcf,'Toolbar','none','Menu','none');  % remove toolbar and menu
            set(gcf,'Name','EEG (blue) and artifacts removed (red)','NumberTitle', 'Off')  % change figure name
            finish_figure(gcf)
        end
    end
    
    % Interpolate bad channels only after ASR: interpolated channels lower the
    % data rank, which hurts ASR's PCA (ICA below is run at the effective rank)
    if EEG.nbchan>10
        EEG = pop_interp(EEG, params.orichanlocs, 'spherical'); % interpolate
        EEG.etc.clean_channel_mask(1:EEG.nbchan) = true;
        if isfield(params,'ref_applied') && strcmp(params.ref_applied,'rest')
            % the original channel info brings back the original reference label
            [EEG.chanlocs.ref] = deal('rest');
            EEG.ref = 'rest';
        end
        if isfield(params,'ica_source') && ~isempty(params.ica_source)
            params.ica_source = pop_interp(params.ica_source, params.orichanlocs, 'spherical');
        end
    else
        warndlg('Cannot interpolate bad EEG channels reliably with less than 10 channels')
        warning('Cannot interpolate bad EEG channels reliably with less than 10 channels')
    end

    % With GEDAI, ICA only for the IC measures of the HEP (nothing removed)
    doICs = strcmp(params.analysis,'hep') && isfield(params,'hep_level') && any(strcmp(params.hep_level,{'ics' 'both'}));
    if useGEDAI && ~doICs
        params.removed_eeg_components = [];
        params.clean_eeg_step = 2;
        return
    end

    % Data the ICA is fitted on: with a high-pass below 1 Hz, a 1-Hz
    % high-passed copy (slow drifts degrade ICA and ICLabel)
    EEGica = EEG;
    if highpass < 1
        if isfield(params,'ica_source') && ~isempty(params.ica_source)
            EEGica = params.ica_source;          % epoched copy from RUN_HEP
        elseif EEG.trials == 1
            EEGica = pop_eegfiltnew(EEG,'locutoff',1,'minphase',causalfilt);
        end
        if size(EEGica.data,1) ~= size(EEG.data,1) || size(EEGica.data,3) ~= size(EEG.data,3)
            warning('The 1-Hz copy for ICA does not match the data: ICA is fitted on the data themselves.')
            EEGica = EEG;
        else
            fprintf('ICA fitted on a 1-Hz high-passed copy of the data. \n')
        end
    end

    % Run ICA at the effective data rank to avoid ghost ICs (Kim et al. 2023).
    % Method 3 uses lrate = 1e-5 and maxsteps = 2000 for reproducible results.
    dataRank = sum(eig(cov(double(EEGica.data(:,:)'))) > 1E-7);
    if icamethod == 1
        EEGica = pop_runica(EEGica,'icatype','picard','maxiter',400,'mode','standard', 'pca',dataRank);
    elseif icamethod == 2
        EEGica = pop_runica(EEGica,'icatype','runica','extended',1,'pca',dataRank);
    elseif icamethod == 3
        EEGica = pop_runica(EEGica,'icatype','runica','extended',1, ...
            'pca',dataRank,'lrate',1e-5,'maxsteps',2000);
    end
    EEG.icaweights = EEGica.icaweights;
    EEG.icasphere = EEGica.icasphere;
    EEG.icawinv = EEGica.icawinv;
    EEG.icachansind = EEGica.icachansind;
    EEG.icaact = [];
    EEG = eeg_checkset(EEG);
    params = rmfield_if(params, 'ica_source');

    % With GEDAI (IC measures): classify the components, remove none
    if useGEDAI
        EEG = pop_iclabel(EEG,'default');
        params.removed_eeg_components = [];
        params.clean_eeg_step = 2;
        return
    end

    % Classify and flag bad components with ICLabel. pop_icflag rows:
    % brain, muscle, eye, heart, line noise, channel noise, other
    EEG = pop_iclabel(EEG,'default');
    if contains(params.analysis, 'rm_heart')
        % Keep heart components: remove_heartcomp removes them next with the
        % user's confidence threshold
        EEG = pop_icflag(EEG,[NaN NaN; .99 1; .9 1; NaN NaN; .99 1; .99 1; NaN NaN]);
    elseif contains(params.analysis, 'hep')
        % HEP: heart (cardiac field artifact) components are removed too
        % with params.heart_removal = 'ica' (default); params.conf_thresh:
        % minimum ICLabel heart probability (default .75)
        conf_thresh = .75;
        if isfield(params,'conf_thresh') && ~isempty(params.conf_thresh), conf_thresh = params.conf_thresh; end
        if isfield(params,'heart_removal') && ~strcmpi(params.heart_removal,'ica'), conf_thresh = NaN; end
        EEG = pop_icflag(EEG,[NaN NaN; .99 1; .9 1; conf_thresh 1; .99 1; .99 1; NaN NaN]);
    else
        % Features and coherence
        EEG = pop_icflag(EEG,[NaN NaN; .95 1; .95 1; .99 1; .99 1; .99 1; NaN NaN]);
    end
    badComp = find(EEG.reject.gcompreject);
    EEG = eeg_checkset(EEG);

    % Store in params if users want that information
    params.removed_eeg_components = badComp;

    % Plot the first 24 independent components (flagged ones marked)
    if params.vis_cleaning
        nComps = size(EEG.icaweights,1);
        if ~isempty(nComps) && nComps>0
            if nComps >= 24
                pop_selectcomps(EEG,1:24);
                set(gcf,'Toolbar','none','Menu','none','Name','Independent components','NumberTitle','Off');  % remove toolbar and menu, set name
            else
                pop_selectcomps(EEG,1:nComps);
                set(gcf,'Toolbar','none','Menu','none','Name','Independent components','NumberTitle','Off');  % remove toolbar and menu, set name
                
            end
            colormap("parula"); finish_figure(gcf)
        else
            warndlg("No independent components in your dataset. Somethig went wrong with your ICA decomposition.")
        end
    end
    
    % Remove bad components
    if ~isempty(badComp)
        fprintf('Removing %g bad component(s). \n', length(badComp));
        EEG = pop_subcomp(EEG, badComp, 0);
    end
    
    % update tracker
    params.clean_eeg_step = 2;

end


function s = rmfield_if(s, f)
% Remove field f from struct s if present
if isfield(s, f), s = rmfield(s, f); end


function EEG = rest_ref(EEG)
% Reference to infinity (REST plugin). The plugin's 'files' folder (dipole
% file, leadfield function) is put on the path here: EEGLAB started with
% EEGLAB does not always add the plugin's subfolders, and REST then fails.
restDir = fileparts(which('ref_infinity'));
if ~isempty(restDir) && exist(fullfile(restDir,'files'),'dir')
    addpath(fullfile(restDir,'files'));
end
EEG = ref_infinity(EEG);
% channel-level label too: eeg_checkset copies chanlocs(1).ref into EEG.ref
% (GEDAI sets it to 'average')
if isfield(EEG.chanlocs,'ref')
    [EEG.chanlocs.ref] = deal('rest');
end
