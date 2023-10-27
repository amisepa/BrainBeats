% BRAINBEATS_PROCESS - Process single EEGLAB files containg EEG and cardiovascular (ECG or PPG) signals.
% 
% Usage:
%    [EEG, com] = brainbeats_process(EEG, 'key', 'val')
% 
% Inputs:
%  'analysis'       - 'hep' (heartbeat-evoked potentials) | 'features' (extract EEG and HRV features) |
%                       'rm_heart' | extract heart components from EEG signals.
%  'heart_signal'   - [ppg'|'ecg'] define if cardiovascular signal to process is PPG or ECG.
%  'heart_channels' - [cell array of character] name(s) of the heart channel to process
%  'rr_correct'     - [see below] correction method for RR artifacts. Interpolation algorithms:
%                       'pchip'     (default; shape-preserving piecewise cubic), 
%                       'linear' 
%                       'cubic'     (falls back to 'spline' interpolation for irregularly-spaced data)
%                       'nearest'   (nearest neighbor)
%                       'next'      (next neighbor)
%                       'previous'  (previous neighbor)
%                       'spline'    (piecewise cubic spline) 
%                       'makima' (modified Akima cubic interpolation)
%                    Or remove them instead with: 'remove'.
%  'clean_eeg'      - [0|1] filter EEG data (bandpass zero-phase FIR 1-45 Hz), 
%                       re-reference data do infinity using REST, remove
%                       bad channels and interpolate them. For HEP, bad
%                       trials are removed, whereas for Features, artifacts
%                       are removed with Artifact Subspace Reconstruction
%                       (ASR). Source separation (ICA) is performed taking
%                       into account effective data rank (see Kim et al.
%                       2023), before classifying and removing bad
%                       components with ICLabel (muscle and heart with 95%
%                       confidence, eye with 90% confidence).
%  'parpool'        - [0|1] use paralell toolbox. Default is 0.
%  'rm_heart'       - [0|1] remove heart channel (1, default) after processing it,
%                     or not (0). 
%  'hrv_features'   - [cell array of characters] HRV features to compute. See
%                     GET_HRV_FEATURES for more information. Choices are
%                     'time' (time-domain measures), 'frequency' (frequency
%                     domain measures, and 'nonlinear' (nonlinear domain measures).
%  'hrv_spec'       - ['LombScargle'|'pwelch'|'fft'] method to compute the 
%                     HRV spectrum. Default is 'LombScargle'. pwelch and
%                     fft implement resampling. 
%  'eeg_features'   - [cell array of string] EEG features to compute. See
%                     GET_EEG_FEATURES for more information. Choices are
%                     'time' (time domain), 'frequency' (frequency domain), 
%                     and 'nonlinear' (nonlinear domain).
%  'norm'           - [0|1] normalize HRV and EEG spectra (Features mode).
%                     For HRV, applied during Lomb-Scargle periodogram
%                     estimation by scaling the total power with variance
%                     in the time series, and in a 2nd step dy dividing
%                     each band power by total power to provide more
%                     intuitive measure of the relative contribution of
%                     each frequency component to overall power. 
%  'gpu'            - [0|1] use GPU. Default is 0.
%  'vis'            - [0|1] set vizualization to on (1) or off (0). Default is on.
%  'save'           - [0|1] save results into a MATLAB file (1) or not (0). Default is on.
%
% Outputs:
%   EEG      - modified EEGLAB dataset. A EEG.brainbeats field is created
%              containing all the measures computed for the dataset.
%   features - Features computed (same as field EEG.brainbeats)
%
% Other options (not documented yet)
% Method 1: Hearbteat evoked potentials (HEP) and oscillations (HEO).
% Method 3: Remove heart components from EEG signals.
% 
% Copyright (C) - Cedric Cannard, 2023

function [EEG, com] = brainbeats_process(EEG, varargin)

Features = [];
com = '';

% Basic checks
if ~exist('EEG','var')
    error('This plugin requires that your data (containing both EEG and ECG/PPG signals) are already loaded into EEGLAB.')
end
if nargin < 1
    help pop_BrainBeats; return;
end
if isempty(EEG) || isempty(EEG.data)
    error('Empty EEG dataset.');
end
if isempty(EEG.chanlocs(1).labels)
    error('No channel labels.');
end
if isempty(EEG.ref)
    warning('EEG data not referenced! Referencing is highly recommended');
end

pop_editoptions('option_single', 0); % ensure double precision

% Add path to subfolders
mainpath = fileparts(which('eegplugin_BrainBeats.m'));
addpath(fullfile(mainpath, 'functions'));
% addpath(fullfile(mainpath, 'functions', 'restingIAF'));
% outPath = fullfile(mainpath, 'sample_data'); %FIXME: ASK USER FOR OUTPUT DIR

% Install dependency plugins
fprintf('Installing necessary EEGLAB plugins... \n')
if ~exist('clean_asr','file')
    plugin_askinstall('clean_asr','clean_asr', 0);
end
if ~exist('picard','file')
    plugin_askinstall('picard', 'picard', 0);
end
if ~exist('iclabel','file')
    plugin_askinstall('iclabel', 'iclabel', 0);
end


%%%%%%%%%%%%%%%%%%%% Main parameters %%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    params = getparams_gui(EEG);                % GUI
elseif nargin > 1
    params = getparams_command(varargin{:});    % Command line
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if GUI was aborted (FIXME: should not send this error)
if ~isfield(params, 'heart_channels')
    disp('Aborted'); return
end

% Check if data format is compatible with chosen analysis and select analysis
if isfield(params,'analysis')
    switch params.analysis
        case 'continuous'
            if length(size(EEG.data)) ~= 2
                error("You selected feature-based analysis but your data are not continuous.")
            end
        case 'epoched'
            if length(size(EEG.data)) ~= 3
                error("You selected HEP analysis but your data are not epoched.")
            end
    end
else
    % Select analysis based on data format if not defined
    if length(size(EEG.data)) == 2
        params.analysis = 'continuous';
        disp("Analysis not defined. Continuous data detected: selecting 'feature-based mode' by default")
    elseif length(size(EEG.data)) == 3
        params.analysis = 'epoched';
        disp("Analysis not defined. Epoched data detected: selecting 'heart-beat evoked potential (HEP) mode' by default")
    else
        error("You did not define the analysis to run, and your data format was not recognized. " + ...
            "Should be 'continuous' or 'epoched', and something may be wrong with your data format ")
    end
end

% Check if heart channels are in file (for command line mode)
if contains({EEG.chanlocs.labels}, params.heart_channels) == 0
    error('The heart channel names you inputted cannot be found in the current dataset.')
else
    % make sure it's a cell if only one ECG channel is provided
    if ~iscell(params.heart_channels)
        params.heart_channels = {params.heart_channels};
    end
end

% Check heart signal type
if ~contains(params.heart_signal, {'ecg' 'ppg'})
    error('Heart signal should be either ECG or PPG')
end

% Check for channel locations for visualization
if params.vis && params.eeg
    if ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(1).X)
        error("Electrode location coordinates must be loaded for visualizing outputs.");
    end
end

% Parallel computing
if ~isfield(params,'parpool') 
    params.parpool = false;
end

% GPU computing
if ~isfield(params,'gpu') 
    params.gpu = false;
end

% Save outputs?
if ~isfield(params,'save') 
    params.save = true;
end

%%%%%%%%%%%%% PREP EEG DATA %%%%%%%%%%%%%
EEG.data = double(EEG.data);  % ensure double precision
params.fs = EEG.srate;

%%%%% MODE 1: remove heart components from EEG signals with IClabel %%%%%
if strcmp(params.analysis,'rm_heart')
    EEG = remove_heartcomp(EEG, params);
end

%%%%% MODE 2 & 3: RR, SQI, and NN %%%%%
if contains(params.analysis, {'features' 'hep'})
    
    % Keep only EEG and ECG channels in separate structures
    % EEG = pop_select(EEG, 'chantype',{'ECG','EEG'});  
    CARDIO = pop_select(EEG,'channel',params.heart_channels); % export ECG data in separate structure
    EEG = pop_select(EEG,'nochannel',params.heart_channels); 

    % Get RR and NN intervals from ECG/PPG signals
    % note: when several electrodes are provided, use the elec with best 
    % quality for subsequent analysis. Structures to avoid issues when 
    % intervals have different across electrodes.
    signal = CARDIO.data;
    nElec = size(signal,1);
    for iElec = 1:nElec
        elec = sprintf('elec%g',iElec);
        fprintf('Detecting R peaks from ECG time series: electrode %g...\n', iElec)
        [RR.(elec), RR_t.(elec), Rpeaks.(elec), sig(iElec,:), sig_t(iElec,:), HR] = get_RR(signal(iElec,:)', params.fs, params.heart_signal);

        % SQI
        SQIthresh = .9; % minimum SQI recommended by Vest et al. (2019)
        if strcmpi(params.heart_signal, 'ecg')
            [sqi(iElec,:), sqi_times(iElec,:)] = get_sqi(Rpeaks.(elec), signal(iElec,:), params.fs);
            SQI(iElec,:) = sum(sqi(iElec,:) < SQIthresh) / length(sqi(iElec,:));  
        % else      
        %     [sqi, ~, annot]= get_sqi_ppg(beats,ppg',params.fs,30); %FIXME: doesn't work
        %     sqi = [beats'./params.fs, sqi'./100];
        end

        % Correct RR artifacts (e.g., arrhytmia, ectopy, noise) to obtain the NN series
        % FIXME: does not take SQI into account
        % tmp_t = RR_t.(elec); 
        % if strcmpi(params.heart_signal, 'ecg')
        % rr_t(1) = [];   % remove 1st heartbeat
        % end
        vis = false;    % more detailed visualization of RR artifacts
        [NN.(elec), NN_t.(elec), flagged.(elec)] = clean_rr(RR_t.(elec), RR.(elec), params, vis);
        flaggedRatio.(elec) = sum(flagged.(elec)) / length(flagged.(elec));

    end

    % Keep only ECG data of electrode with the lowest number of RR
    % artifacts (FIXME: add SQI)
    % [~, best_elec] = min(SQI);
    % sqi = [sqi_times(best_elec,:); sqi(best_elec,:)];
    % SQI = SQI(best_elec);
    [~,best_elec] = min(struct2array(flaggedRatio));
    elec = sprintf('elec%g',best_elec);
    maxThresh = .2;         % max portion of artifacts (.2 default from Vest et al. 2019)
    % if SQI > maxThresh    % more than 20% of RR series is bad
    %      warning on
    %     warning("%g%% of the RR series on your best ECG electrode has a signal quality index (SQI) below minimum recommendations (max 20%% below SQI = .9; see Vest et al., 2019)!",round(SQI,2));
    %     error("Signal quality is too low: aborting! You could inspect the data in EEGLAB > Plot > Channel data (Scroll) and try to remove large artifacts first.");
    % else
    %     fprintf( "Keeping only the heart electrode with the best signal quality index (SQI): %g%% of the RR series is outside of the recommended threshold. \n", SQI )
    % end
    flaggedRatio = flaggedRatio.(elec);
    flagged = flagged.(elec);
    if sum(flagged) > 0
        if contains(params.rr_correct,'remove')
            fprintf('%g heart beats were flagged as artifacts and removed. \n', sum(flagged));
        else
            fprintf('%g heart beats were flagged as artifacts and interpolated. \n', sum(flagged));
        end
    end
    if flaggedRatio > maxThresh % more than 20% of RR series is bad
        warning("%g%% of the RR series on your best electrode are artifacts. Maximum recommendations 20%%! Aborting...", round(flaggedRatio*100,2));
    else
        fprintf( "Keeping only the heart electrode with the best signal quality index (SQI): %g%% of the RR series is outside of the recommended threshold. \n", round(flaggedRatio,2) )
    end
    sig_t = sig_t(best_elec,:);
    sig = sig(best_elec,:);
    RR = RR.(elec);
    RR_t = RR_t.(elec);
    RR_t(1) = [];       % always ignore 1st hearbeat
    Rpeaks = Rpeaks.(elec);
    Rpeaks(1) = [];     % always ignore 1st hearbeat
    NN_t = NN_t.(elec);
    NN = NN.(elec);

    % Plot filtered ECG and RR series of best electrode and interpolated
    % RR artifacts (if any)
    if params.vis
        plot_NN(sig_t, sig, RR_t, RR, Rpeaks, NN_t, NN, flagged)
    end
    
    % Preprocessing outputs
    Features.HRV.ECG_filtered = sig;
    Features.HRV.ECG_times = sig_t;
    if exist('SQI','var')
        Features.HRV.SQI = SQI;
    end
    Features.HRV.RR = RR;
    Features.HRV.RR_times = RR_t;
    Features.HRV.HR = HR;
    Features.HRV.NN = NN;
    Features.HRV.NN_times = NN_t;
    Features.HRV.flagged_heartbeats = flagged;

    % Remove ECG data from EEG data
    EEG = pop_select(EEG,'nochannel',params.heart_channels); % FIXME: remove all non-EEG channels instead

    % Filter, re-reference, remove bad channels
    if params.clean_eeg    
        params.clean_eeg_step = 0;
        [EEG, params] = clean_eeg(EEG, params);
    end

    %%%%% MODE 2: Heartbeat-evoked potentials (HEP) %%%%%
    if strcmp(params.analysis,'hep')
        EEG = run_HEP(EEG, params, Rpeaks);
    end 
    
    %%%%% MODE 3: HRV features %%%%%
    if strcmp(params.analysis,'features') && params.hrv

            % File length (we take the whole series for now to allow ULF and VLF as much as possible)
            file_length = floor(EEG.xmax)-1;
            if file_length < 300
                warning('File length is shorter than 5 minutes! The minimum recommended is 300 s for estimating reliable HRV metrics.')
            end
            params.file_length = file_length;

            % Extract HRV measures
            HRV = get_hrv_features(NN, NN_t, params);

            % Final output with everything
            Features.HRV = HRV;
    end

    %%%%% MODE 3: EEG features %%%%%
    if strcmp(params.analysis,'features') && params.eeg

        % Clean EEG artifacts with ASR
        if params.clean_eeg
            [EEG, params] = clean_eeg(EEG, params);
        end

        % Extract EEG features and store in EEGLAB structure
        if params.eeg
            params.chanlocs = EEG.chanlocs;
            EEG.features = get_eeg_features(EEG.data, params);
        end

        % Final output with everything to export separately as .mat file
        Features.EEG = EEG.features;

    end

end

%%%%%% PLOT & SAVE FEATURES %%%%%%%
if strcmp(params.analysis,'features')
    
    % Save in same repo as loaded file (FIXME: ASK USER FOR OUTPUT DIR)
    if params.save
        outputPath = fullfile(EEG.filepath, sprintf('%s_features.mat', EEG.filename(1:end-4)));
        fprintf("Saving features in EEG.features and exporting them here: %s \n", outputPath);
        save(outputPath,'Features');
    end

    % Plot features
    if params.vis
        plot_features(Features,params)
    end
end

% Shut down parallel pool
% if params.parpool
%     delete(gcp('nocreate'));
% end


%%%%%%%%%%%%%%%%%%%%%%%% eegh %%%%%%%%%%%%%%%%%%%%%%%%
if contains(params.analysis,'hep')
    com = char(sprintf("EEG = brainbeats_process(EEG,'analysis','%s','heart_signal','%s','heart_channels','%s','rr_correct','%s','clean_eeg',%g,'parpool',%g,'gpu',%g,'vis',%g,'save',%g);",...
        params.analysis,params.heart_signal,['{' sprintf(' ''%s'' ', params.heart_channels{:}) '}'],params.rr_correct,params.clean_eeg,params.parpool,params.gpu,params.vis,params.save));
elseif contains(params.analysis,'features') 
    com = char(sprintf("[~, Features] = brainbeats_process(EEG,'analysis','%s','heart_signal','%s','heart_channels','%s','rr_correct','%s','clean_eeg',%g,'parpool',%g,'gpu',%g,'vis',%g,'save',%g);",...
        params.analysis,params.heart_signal,['{' sprintf(' ''%s'' ', params.heart_channels{:}) '}'],params.rr_correct,params.clean_eeg,params.parpool,params.gpu,params.vis,params.save));
elseif contains(params.analysis,'rm_heart') 
    com = char(sprintf("EEG = brainbeats_process(EEG,'analysis','%s','heart_signal','%s','heart_channels','%s','clean_eeg',%g,'vis',%g,'save',%g);",...
        params.analysis,params.heart_signal,['{' sprintf(' ''%s'' ', params.heart_channels{:}) '}'],params.clean_eeg,params.vis,params.save));
end
com = char(com);

%%%%%%%%%%%%%%%%%%%%%%%% References %%%%%%%%%%%%%%%%%%%%%%%%
if params.hrv
    fprintf("For cardivosacular signal processing and HRV metrics, please cite: \n ");
    fprintf("   - Vest et al. (2018). An open source benchmarked toolbox for cardiovascular waveform and interval analysis. Physiol Meas. \n");
    fprintf("   - Shaffer & Ginsberg (2017). An overview of heart rate variability metrics and norms. Frontiers in public health. \n");
end
if strcmp(params.analysis,'features')
%     if params.eeg && params.eeg_frequency
%         fprintf("For the IAF feature, please cite: \n %s \n", "  - Corcoran et al. (2018). Toward a reliable, automated method of individual alpha frequency (IAF) quantification. Psychophysiology. ")
%         fprintf("For the alpha asymmetry feature, please cite: \n %s \n", "  - Smith et al. (2017). Assessing and conceptualizing frontal EEG asymmetry: An updated primer on recording, processing, analyzing, and interpreting frontal alpha asymmetry. International Journal of Psychophysiology.")
%         fprintf("For the EEG coherence measures, please cite: \n");
%         fprintf("   - Cao et al. (2016). Resting-state EEG power and coherence vary between migraine phases. The journal of headache and pain. \n");
%         fprintf("   - Pascual-Marqui et al. (2014). Assessing direct paths of intracortical causal information flow of oscillatory activity with the isolated effective coherence (iCoh). Frontiers in human neuroscience. \n");
%     end
    if params.hrv_nonlinear || params.eeg_nonlinear
        fprintf("For the entropy measures, please cite: \n");
        fprintf("   - Azami & Escudero (2016). Refined Multiscale Fuzzy Entropy based on standard deviation for biomedical signal analysis. Medical & Biological Engineering & Computing. \n");
%         fprintf("   - Costa, Goldberger, Peng (2002). Multiscale entropy analysis of complex physiologic time series. Phys Rev Lett. \n ")
%         fprintf("   - Kosciessa et al. (2020). Standard multiscale entropy reflects neural dynamics at mismatched temporal scales: What's signal irregularity got to do with it? Plos Comput Biol. \n ")
    end
end
% if strcmp(params.analysis,'hep')
%     fprintf("For HEP methods, please cite: \n");
%     fprintf("   - Candia-Rivera et al. (2021). The role of electroencephalography electrical reference in the assessment of functional brain–heart interplay: From methodology to user guidelines. Neuroscience Methods. \n");
%     fprintf("   - Park & Blanke (2019). Heartbeat-evoked cortical responses: Underlying mechanisms, functional roles, and methodological considerations. NeuroImage. \n");
% end
fprintf('\n\n')
fprintf("Thank you for using the BrainBeats toolbox! \nPlease cite: \nCannard, Wahbeh, & Delorme (2023). BrainBeats: an open-source EEGLAB plugin to jointly analyze EEG and cardiovascular (ECG/PPG) signals. bioRxiv, 2023-06 \n")

disp('Done!'); gong
