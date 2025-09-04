%% BrainBeats clean_eeg 
% 
% 1) Applies a lowpass at 1 Hz and highpass at 40 Hz are applied.
%   For HEP, a nonlinear minimum-phase FIR filter is used to preserve 
%   causality (especially important for pre-heartbeat analysis, whereas a 
%   zero-phase noncausal FIR filter is used for continuous data.
% 2) If data were not already referenced and have at least 30 channels, 
%   they re-referenced to average. 
% 3) bad EEG channels are identified and reomved using the clean_flatlines, 
%   and clean_channels algorithms from K. Kothe (clean_artifacts). 
%   Default parameters: 
%       - correlation threshold = .65; 
%       - window length = 5 s to capture low-frequency artifacts and increase 
%           speed; 
%       - line noise threshold = 15; 
%       - maximum portion of channel to be considred a bad channel = 33%; 
%       - 85% of available RAM to increase speed.
%       - # of ransac samples = 500 (more computation but more reliable and 
%           replicable)
% 4) Bad channels are interpolated using spherical splines. 
% 5a) For continuous data, large artifacts are removed using Artifact 
%   subspace reconstruction (ASR). SD threshold = 30 by default.
% 5b) For epoched data (HEP), bad epochs are detected and removed using
%   custom amplitude and SNR metrics. Default method = 'grubbs'.
% 6) The Infomax algorithm is used by default, implementing PCA-dimension 
%   reduction to effectiv edata rank to avoid ghost ICs (see Kim et al. 2023). 
%  Set 'icamethod' to 1 to use the PICARD algorithm (much faster but
%  requires installation of the plugin, which does't always work automatically). 
%   Set 'icamethod' to 3 to use the modified Infomax ('lrate' = 1e-5 and 
%   'maxsteps' = 2000) to increase convergence/replicability (but takes
%   much longer). 
% 7) classify ICs with ICLabel: 
%       - eye with 90% confidence
%       - muscle, heart, line noise, channel noise with with 99% confidence
%           NOTE: heart is not removed for 'hep' and 'rm_heart' methods 
%               since we want to preserve that activity.
% 
% Copyright (C), BrainBeats 2023, Cedric Cannard

function [EEG, params] = clean_eeg(EEG, params)

% General parameters
if isfield(params,'ref')
    reref = params.ref;
else
    reref = 'average'; % 'average' (default), 'infinity', 'csd', 'off'
    params.ref = 'average';    % for user
end
if isfield(params,'highpass')
    highpass = params.highpass;
else
    highpass = 1; % default = 1 Hz
    params.highpass = 1;    % for user
end
if isfield(params,'lowpass')
    lowpass = params.lowpass;
else
    lowpass = 40;   % default = 40 hz to safely remove power line at 50/60 Hz 
                    % with a low filter order (less ripples and faster)
    params.lowpass = 40;    % for user
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
win_length = 2;     % window length to scan channels (default = 2 s)
line_thresh = 15;    % line noise threshold to remove bad channels (default = 15)
nSamp = 200;        % number of ransac samples (default = 200; ~50-500 range; higher is longer but more accurate and replicable)

% HEP parameters to remove bad epochs
if isfield(params,'detectMethod')
    detectMethod = params.detectMethod;
else
    detectMethod = 'grubbs';   % 'median' (more agressive), 'grubbs' (moderate; default), 'mean' (more lax)
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
    icamethod = 2;  % 1 = fast ICA (picard), 2 = Infomax, 3 replicable Infomax (longest but replicable)
end

% Filter, re-reference, and remove bad channels
if params.clean_eeg_step == 0
    
    % Highpass filter to remove slow frequency drifts set by user
    EEG = pop_eegfiltnew(EEG,'locutoff',highpass,'minphase',causalfilt);    

    % Lowpass filter set by user
    EEG = pop_eegfiltnew(EEG,'hicutoff',lowpass,'minphase',causalfilt);  
    
    % Reference to average or infinity/REST
    % Candia-Rivera, Catrambone, & Valenza (2021). The role of EEG reference 
    % in the assessment of functional brain–heart interplay: From 
    % methodology to user guidelines. Journal of Neuroscience Methods.
    if ~strcmp(reref,'off')
        if EEG.nbchan < 30
            warndlg('Cannot reference these EEG data to infinity or average or Surface Laplacian, not validated with less than 30 channels.')
            warning('Cannot reference these EEG data to infinity or average or Surface Laplacian, not validated with less than 30 channels.')
        else
            if strcmp(reref,'infinity')
                fprintf('Re-referencing EEG data to infinity. \n')
                try
                    EEG = ref_infinity(EEG);
                catch
                    warning("Re-reference to infinity failed. This may happen on MACs (you likely need XCode installed for compiling the required code). Please submit an issue on Github: https://github.com/amisepa/BrainBeats/issues")
                    warning("Defaulting back to common average reference (CAR).")
                    EEG = apply_car(EEG);  % preserving effective data rank
                end
            elseif strcmp(reref,'average')
                fprintf('Re-referencing EEG data to average. \n')
                % EEG = pop_reref(EEG,[]);
                EEG = apply_car(EEG);  % preserving effective data rank
            elseif strcmp(reref,'csd')
                try
                    disp("Performing reference-free current-source density (CSD) transformation (Surface Laplacian).")
                    EEG = csd_transform(EEG);
                catch
                    warning("Re-reference to infinity failed. This may happen on MACs. Please submit an issue on Github: https://github.com/amisepa/BrainBeats/issues")
                    warning("Defaulting back to common average reference (CAR).")
                    EEG = apply_car(EEG);  % preserving effective data rank
                end
            end
        end
    end
    
    % Remove bad channels
    EEG.etc.clean_channel_mask = true(1,EEG.nbchan);
    oriEEG = EEG;
    % EEG = pop_clean_rawdata(EEG,'FlatlineCriterion',5,'ChannelCriterion',.85, ...
    %     'LineNoiseCriterion',5,'Highpass','off', 'BurstCriterion','off', ...
    %     'WindowCriterion','off','BurstRejection','off','Distance','off');    
    % disp('Scanning EEG channels for flat lines...')
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
        % EEG.etc.clean_channel_mask(1:EEG.nbchan) = true;
        % EEG.etc.clean_channel_mask(badChan) = false;
        try
            vis_artifacts(EEG,oriEEG,'ShowSetname',false); pause(0.01)
        catch
            warning('failed to plot bad channels or artifacts with vis_artifacts(). Setting show_events to off and trying againg')
            try
                vis_artifacts(EEG,oriEEG,'ShowSetname',false,'ShowEvents',false); pause(0.01)
            catch
                warning("vis_artifacts failed to plot the removed channels or artifacts. Please submit an issue on EEGLAB's page: https://github.com/sccn/eeglab/issues")
            end
        end
        try icadefs; set(gcf, 'color', BACKCOLOR); catch; end     % eeglab background color
        set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
        set(gcf,'Name','EEG channels removed','NumberTitle', 'Off')  % change figure name
        % vis_artifacts(EEG,oriEEG,'ChannelSubset',1:EEG.nbchan-length(params.heart_channels));
    else
        disp("No bad channels detected.")
    end
    
    % notch filter if line noise is below lowpass filter
    if isfield(params,'linenoise') && params.linenoise<lowpass
        EEG = pop_eegfiltnew(EEG, 'locutoff',params.linenoise-3, ...
            'hicutoff',params.linenoise+3,'revfilt',1,'filtorder',500);
    end

    % update tracker
    params.clean_eeg_step = 1;

% Remove bad trials for HEP, aritfacts for Features
elseif params.clean_eeg_step == 1
    
    disp('----------------------------------------------')
    fprintf('              Cleaning EEG data \n')
    disp('----------------------------------------------')

    % HEP (remove bad epochs)
    if strcmp(params.analysis, 'hep')
        
        % Detect and remove bad epochs
        badTrials = find_badTrials(EEG, detectMethod, params.vis_cleaning);
        EEG = pop_rejepoch(EEG, badTrials, 0); 
        
        % Run RMS a 2nd time more conservative in case some were missed
        % badTrials = find_badTrials(EEG,'mean', params.vis_cleaning);
        % EEG = pop_rejepoch(EEG, badTrials, 0);
        
        % Store in params if users want that information
        params.removed_eeg_trials = badTrials;
        
    % ASR on continuous data
    elseif contains(params.analysis, {'features' 'rm_heart' 'coherence'})
        
        % Identify artifacts using ASR
        oriEEG = EEG;
        try
            m = memory; maxmem = round(asr_mem*(m.MemAvailableAllArrays/1000000),1);  % use 80% of available memory (in MB)
            cleanEEG = clean_asr(EEG,asr_cutoff,[],[],[],[],[],[],usegpu,false,maxmem);
        catch
            warning("Failed to use high RAM to run ASR faster. Defaulting back to default values (ASR will just be slower).")
            cleanEEG = clean_asr(EEG,asr_cutoff,[],[],[],[],[],[],usegpu,false,[]);
        end
        
        % Mask for vis_artifacts
        mask = sum(abs(EEG.data-cleanEEG.data),1) > 1e-10;
        EEG.etc.clean_sample_mask = true(1, length(mask)); % initialize all samples as clean
        badData = reshape(find(diff([false mask false])), 2, [])';
        badData(:, 2) = badData(:, 2) - 1;
        % exclude very short artifacts < 10 samples
        if ~isempty(badData)  
            smallIntervals = diff(badData')' < 10;
            badData(smallIntervals, :) = [];
            for i = 1:size(badData, 1)
                EEG.etc.clean_sample_mask(badData(i, 1):badData(i, 2)) = false;
            end
        end

        % Remove them from data
        EEG = pop_select(EEG,'nopoint',badData);
        % if strcmp(params.analysis,'hep')
        %     CARDIO = pop_select(CARDIO,'nopoint',badData);
        % end
        fprintf('%g %% of data were considered to be artifacts and were removed. \n', (1-EEG.xmax/oriEEG.xmax)*100)
        
        % Store in params if users want that information
        params.removed_eeg_segments = badData;

        % Plot what has been removed
        if params.vis_cleaning
            vis_artifacts(EEG,oriEEG,'ShowSetname',false); pause(0.01)
            try icadefs; set(gcf, 'color', BACKCOLOR); catch; end     % eeglab background color
            set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
            set(gcf,'Name','EEG (blue) and artifacts removed (red)','NumberTitle', 'Off')  % change figure name
        end
    end
    
    % Interpolate bad channels (after ASR as low data rank can cause bad
    % performance with PCA used in ASR. But not a problem for ICA as we
    % input the data rank: see below). 
    if EEG.nbchan>10
        EEG = pop_interp(EEG, params.orichanlocs, 'spherical'); % interpolate
        EEG.etc.clean_channel_mask(1:EEG.nbchan) = true;
    else
        warndlg('Cannot interpolate bad EEG channels reliably with less than 10 channels')
        warning('Cannot interpolate bad EEG channels reliably with less than 10 channels')
    end

    % Run ICA at effective data rank to control for ghost ICs (Kim et al. 2023). 
    % Use Picard algorithm by default to increase speed. 
    % use lrate=1e-5 and maxsteps=2000 to obtain reproducible ICA results
    dataRank = sum(eig(cov(double(EEG.data(:,:)'))) > 1E-7);
    if icamethod == 1
        EEG = pop_runica(EEG,'icatype','picard','maxiter',500,'mode','standard', ...
            'pca',dataRank);
    elseif icamethod == 2
        EEG = pop_runica(EEG,'icatype','runica','extended',1,'pca',dataRank);
    elseif icamethod == 3 
        EEG = pop_runica(EEG,'icatype','runica','extended',1, ...
            'pca',dataRank,'lrate',1e-5,'maxsteps',2000);
    end
    
    % Classify and remove bad components with IClabel
    EEG = pop_iclabel(EEG,'default');
    if contains(params.analysis, {'rm_heart' 'hep'}) 
        % % Do not remove heart components here (we clean the EEG only from 
        % % other artifacts)
        conf_thresh = .75;  % confidence threshold for removing CFA
        EEG = pop_icflag(EEG,[NaN NaN; .99 1; .9 1; conf_thresh 1; .99 1; .99 1; NaN NaN]);
        % EEG = pop_icflag(EEG,[NaN NaN; .95 1; .9 1; NaN NaN; .99 1; .99 1; NaN NaN]);
    else
        % Remove components: brain,  muscle, eye, heart, line noise, channel noise, other
        EEG = pop_icflag(EEG,[NaN NaN; .95 1; .95 1; .99 1; .99 1; .99 1; NaN NaN]);
    end
    badComp = find(EEG.reject.gcompreject);
    EEG = eeg_checkset(EEG);

    % Store in params if users want that information
    params.removed_eeg_components = badComp;

    % Visualize indepent components tagged as bad
    if params.vis_cleaning
        nComps = size(EEG.icaweights,1);
        if ~isempty(nComps) && nComps>0
            if nComps >= 24
                pop_selectcomps(EEG,1:24); pause(0.01)
                set(gcf,'Toolbar','none','Menu','none','Name','Independent components','NumberTitle','Off');  % remove toolbobar and menu and name
            else
                pop_selectcomps(EEG,1:nComps); pause(0.01)
                set(gcf,'Toolbar','none','Menu','none','Name','Independent components','NumberTitle','Off');  % remove toolbobar and menu and name
                
            end
            colormap("parula"); pause(0.01)
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
