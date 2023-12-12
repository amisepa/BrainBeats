%% BrainBeats clean_eeg
% When users choose to preprocess their EEG data with BrainBeats:
% 1) If EEG data have not been filtered, a lowpass at 1 Hz and highoass at 
%   45 Hz are applied (linear zero-phase FIR filters). 
% 2) If data were not already referenced and have at least 30 channels, 
%   they re-referenced to infinity (REST method). 
% 3) bad EEG channels are identified and reomved using the clean_flatlines, 
%   and clean_channels algorithms from K. Kothe (clean_artifacts). 
%   Default parameters: 
%       - correlation threshold = .65; 
%       - window length = 5 s to capture low-frequency artifacts; 
%       - line noise threshold = 5; 
%       - maximum portion of channel to be considred a bad channel = 33%; 
%       - 85% of available RAM to increase speed.
%       - # of ransac samples = 500 (more computation but more accurate and 
%       replicable)
% 4) Bad channels are interpolated using spherical splines. 
% 5a) For continuous data, large artifacts are removed using Artifact 
%   subspace reconstruction (ASR). SD threshold = 30 by default.
% 5b) For epoched data (HEP), bad epochs are detected and removed using
%   custom amplitude and SNR metrics. Default method = 'grubbs'.
% 6) The preconditioned ICA algorithm is used to to perform fast but reliable 
%   blind source spearation, implementing PCA-dimension reduction to control 
%   for effective data rank and avoid ghost IC artifacts (see Kim et al. 2023). 
%   If users do not wish to use PICARD, the infomax algorithm is used with 
%   'lrate' = 1e-5 and 'maxsteps' = 2000. It takes much longer, but allows
%   replication (and the order of the IC maps is always the same). 
% 7) ICLabel tags artifactual ICs with: 
%       - muscle (95% confidence);
%       - eye (90% confidence)
%       - heart (99% confidence) --> NOT removed for 'hep' and 'rm_heart' 
%           modes since we want to preserve that activity.
%       - line noise (99% confidence)
%       - channel noise (99% confidence)
%   Flagged components are automatically extracted from the EEG signals. 
% 
% Copyright (C), BrainBeats 2023, Cedric Cannard

function [EEG, params] = clean_eeg(EEG, params)

% General parameters
if isfield(params,'reref')
    reref = params.reref;
else
    reref = 'infinity'; % default = 'infinity'
end
if isfield(params,'highpass')
    highpass = params.highpass;
else
    highpass = 1; % default
end
if isfield(params,'lowpass')
    lowpass = params.lowpass;
else
    lowpass = 40; % default
end
if isfield(params,'filttype')
    if strcmp(params.filttype,'causal')
        causalfilt = true;  % causal (nonlinear) minimum-phase filter
    else 
        causalfilt = false; % zero-phase (linear) noncausal filter
    end
else
    % default filter depending on analysis
    if strcmp(params.analysis,'hep')
        causalfilt = true; % causal minimum-phase filter should be used for HEP to avoid group delays or other filter artifacts
    else
        causalfilt = false; % zero-phase noncausal filter for continuous data
    end
end
if isfield(params,'gpu')
    usegpu = params.gpu;
else
    usegpu = false; % default
end

% Channel removal parameters
if isfield(params,'flatline')
    flatline = params.flatline;
else
    flatline = 5; % max flat segment to remove channel (default = 5 s)
end
if isfield(params,'corrThresh')
    corrThresh = params.corrThresh;
else
    corrThresh = .65;   % correlation threshold to be considered bad (default = .65)
end
if isfield(params,'maxBad')
    maxBad = params.maxBad;
else
    maxBad = .33;       % max tolerated portion of channel to be bad before removal (default = .33)
end
win_length = 5;     % window length to scan channels (default = 5 s)
line_thresh = 5;    % line noise threshold to remove bad channels (default = 5)
nSamp = 500;        % number of ransac samples (default = 500; higher is longer but more accurate and replicable)

% HEP parameters to remove bad epochs
if isfield(params,'detectMethod')
    detectMethod = params.detectMethod;
else
    detectMethod = 'grubbs';   % 'median' (agressive), 'grubbs' (moderate, default), 'mean' (conservative)
end

% ASR parameters
if isfield(params,'asr_cutoff')
    asr_cutoff = params.asr_cutoff;
else
    asr_cutoff = 30;  % main ASR SD cutoff (lower = more aggressive, higher = more lax)
end
if isfield(params,'asr_mem')
    asr_mem = params.asr_mem;
else
    asr_mem = .8;     % available RAM to use for ASR (.85 = 85% of available RAM)
end

% ICA parameters
if isfield(params,'icamethod')
    icamethod = params.icamethod;
else
    icamethod = 1;  % 1 = fast ICA (picard), 2 = normal Infomax, 3 replicable Infomax (slowest but replicable)
end

% Filter, re-reference, and remove bad channels
if params.clean_eeg_step == 0
    
    % Highpass filter to remove slow frequency drifts set by user
    EEG = pop_eegfiltnew(EEG,'locutoff',highpass,'minphase',causalfilt);    

    % Lowpass filter to remove very high-frequency artifacts (for ASR, lowpass set by user is later below)
    EEG = pop_eegfiltnew(EEG,'hicutoff',100,'minphase',causalfilt);  
    
    % Reference to average or infinity/REST
    % Candia-Rivera, Catrambone, & Valenza (2021). The role of EEG reference 
    % in the assessment of functional brain–heart interplay: From 
    % methodology to user guidelines. Journal of Neuroscience Methods.
    if ~strcmp(reref,'off')
        if EEG.nbchan < 30
            warndlg('Cannot reference these EEG data, not recommended with less than 30 channels.')
            warning('Cannot reference these EEG data, not recommended with less than 30 channels.')
        else
            if strcmp(reref,'infinity')
                fprintf('Re-referencing EEG data to infinity... \n')
                EEG = ref_infinity(EEG);
            elseif strcmp(reref,'average')
                EEG = pop_reref(EEG,[]);
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
        EEG = clean_channels(EEG,corrThresh,line_thresh,win_length,maxBad,nSamp); 
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
    if params.vis_cleaning
        % EEG.etc.clean_channel_mask(1:EEG.nbchan) = true;
        % EEG.etc.clean_channel_mask(badChan) = false;
        vis_artifacts(EEG,oriEEG,'ShowSetname',false);
        set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
        set(gcf,'Name','EEG channels removed','NumberTitle', 'Off')  % change figure name
        % vis_artifacts(EEG,oriEEG,'ChannelSubset',1:EEG.nbchan-length(params.heart_channels));
    end

    % Lowpass filter set by user
    EEG = pop_eegfiltnew(EEG,'hicutoff',lowpass,'minphase',causalfilt);  

    % Interpolate them
    if EEG.nbchan>20
        EEG = pop_interp(EEG, oriEEG.chanlocs, 'spherical'); % interpolate
        EEG.etc.clean_channel_mask(1:EEG.nbchan) = true;
    else
        warning('Cannot interpolate bad EEG channels reliably with less than 30 channels')
    end

    % % Add CARDIO channels back? (requires adding CARDIO in inputs)
    % EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
    % EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
    % for iChan = 1:CARDIO.nbchan
    %     EEG.chanlocs(end+1).labels = params.heart_channels{iChan};
    % end
    % EEG = eeg_checkset(EEG);

    % update tracker
    params.clean_eeg_step = 1;

% Remove bad trials for HEP, aritfacts for Features
elseif params.clean_eeg_step == 1
    
    disp('----------------------------------------------')
    fprintf('              Cleaning EEG data \n')
    disp('----------------------------------------------')

    % HEP
    if strcmp(params.analysis, 'hep')
        
        % Detect and remove bad epochs
        badTrials = find_badTrials(EEG, detectMethod, params.vis_cleaning);
        EEG = pop_rejepoch(EEG, badTrials, 0); 

        % Run RMS a 2nd time more conservative in case some were missed
        % badTrials = find_badTrials(EEG,'mean', params.vis_cleaning);
        % EEG = pop_rejepoch(EEG, badTrials, 0);

        % Store in params if users want that information
        params.removed_eeg_trials = badTrials;

        % Features
    elseif contains(params.analysis, {'features' 'rm_heart'})

        % Identify artifacts using ASR
        oriEEG = EEG;
        m = memory; maxmem = round(asr_mem*(m.MemAvailableAllArrays/1000000),1);  % use 80% of available memory (in MB)
        cleanEEG = clean_asr(EEG,asr_cutoff,[],[],[],[],[],[],usegpu,false,maxmem);
        mask = sum(abs(EEG.data-cleanEEG.data),1) > 1e-10;
        EEG.etc.clean_sample_mask = ~mask;
        badData = reshape(find(diff([false mask false])),2,[])';
        badData(:,2) = badData(:,2)-1;

        % Ignore very small artifacts (<5 samples)
        if ~isempty(badData)
            smallIntervals = diff(badData')' < 5;
            badData(smallIntervals,:) = [];
        end

        % Remove them
        EEG = pop_select(EEG,'nopoint',badData);
        % if strcmp(params.analysis,'hep')
        %     CARDIO = pop_select(CARDIO,'nopoint',badData);
        % end
        fprintf('%g %% of data were considered to be artifacts and were removed. \n', (1-EEG.xmax/oriEEG.xmax)*100)
        
        % Store in params if users want that information
        params.removed_eeg_segments = badData;

        % Plot what has been removed
        if params.vis_cleaning
            vis_artifacts(EEG,oriEEG,'ShowSetname',false);
            set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
            set(gcf,'Name','EEG artifacts channels removed','NumberTitle', 'Off')  % change figure name
        end
    end
    
    % Run ICA at effective data rank to control for ghost ICs (Kim et al. 2023). 
    % Use Picard algorithm when possible to increase speed. 
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
        % Do not remove heart components
        EEG = pop_icflag(EEG,[NaN NaN; .95 1; .9 1; NaN NaN; .99 1; .99 1; NaN NaN]);
    else
        % Remove components: muscle, eye, heart, line noise, channel noise
        EEG = pop_icflag(EEG,[NaN NaN; .95 1; .9 1; .99 1; .99 1; .99 1; NaN NaN]);
    end
    badComp = find(EEG.reject.gcompreject);
    EEG = eeg_checkset(EEG);

    % Store in params if users want that information
    params.removed_eeg_components = badComp;

    % Visualize indepent components tagged as bad
    if params.vis_cleaning
        nComps = size(EEG.icaact,1);
        if nComps >= 24
            pop_selectcomps(EEG,1:24); 
            set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
            set(gcf,'Name','Independent components','NumberTitle','Off')  % name
        else
            pop_selectcomps(EEG,1:nComps)
            set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
            set(gcf,'Name','Independent components','NumberTitle','Off')  % name
        end
        colormap("parula")
    end
    
    % Remove bad components
    if ~isempty(badComp)
        fprintf('Removing %g bad component(s). \n', length(badComp));
        EEG = pop_subcomp(EEG, badComp, 0);
    end
    
    % plot final cleaned data
    if params.vis_cleaning
        if strcmp(params.analysis, 'hep')
            pop_eegplot(EEG,1,1,1);
            set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
            set(gcf,'Name','Final output after preprocessings','NumberTitle','Off')  % name
        % elseif strcmp(params.analysis, 'features')
        %     eegplot(EEG.data,'winlength',15,'srate',EEG.srate,'events', ...
        %         EEG.event,'spacing',50);
        end
    end

    % update tracker
    params.clean_eeg_step = 2;

end 
