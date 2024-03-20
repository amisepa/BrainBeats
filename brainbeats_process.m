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
%                       Or remove them with 'remove' (not recommended).
%  'clean_eeg'      - [0|1] to preprocess EEG data (1) or not (0).
%                       - filter EEG data (bandpass zero-phase FIR 1-45 Hz)
%                       - re-reference data do average (default)
%                       - remove bad channels and interpolate them.
%                       - remove bad epochs (HEP-mode) or bad segments
%                           (Features-mode), whereas for Features
%                       - extract bad components with ICA, taking into
%                           account effective data rank (see Kim et al.2023)
%                           and classifying them with ICLabel (eye with 90%
%                           confidence, others with 99% confidence).
%  'rm_heart'       - [0|1] remove the heart channel (1, default) after
%                       operations are complete (1, default) or not (0).
%  'hrv_features'   - [cell array of characters] HRV features to compute.
%                       Choices are 'time', 'frequency', 'nonlinear'.
%                       See GET_HRV_FEATURES for more information.
%  'eeg_features'   - [cell array of string] EEG features to compute.
%                       Choices are 'frequency', 'nonlinear'.
%                       See GET_EEG_FEATURES for more information.
%  'parpool'        - [0|1] use parallel computing (1) or not (0, default).
%                       Mainly useful for nonlinear EEG features.
%  'gpu'            - [0|1] use GPU computing (1) or not (0).
%                       Mainly useful for cleaning EEG with ASR faster.
%  'vis_cleaning'   - [0|1] vizualize the preprocessing plots (1, default)
%                     or not (0). Strongly recommended.
%  'vis_outputs'    - [0|1] vizualize the outputs (1, default) or not (0).
%                       Include for example grand average HEP/HEO, EEG/HRV
%                       power spectral density, EEG features topographies.
%  'save'           - [0|1] save the resulting processed EEGLAB dataset (.set
%                       file at the same place as loaded file but renamed).
%                       For features-mode, also saves a .mat file containing
%                       the features in structure format. Set to ON (1, default)
%                       or OFF (0).
%
% Outputs:
%   'EEG'           - Processed EEGLAB dataset. An EEG.brainbeats field is
%                   created containing the
%                   parameters used, some preprocessing outputs (e.g. channels
%                   or RR artifacts removed), and the computed EEG & HRV
%                   features computed if using the features-mode.
%
% Copyright (C) - BrainBeats, Cedric Cannard, 2023

function [EEG, com] = brainbeats_process(EEG, varargin)

Features = [];
com = '';

% Basic checks
if ~exist('EEG','var')
    errordlg('The BrainBeats plugin requires that your data (containing both EEG and ECG/PPG signals) are already loaded into EEGLAB.')
end
if nargin < 1
    help brainbeats_process; return;
end

% Add path to subfunctions
mainpath = fileparts(which('eegplugin_BrainBeats.m'));
addpath(fullfile(mainpath, 'functions'));

% basic checks on EEG data
if isempty(EEG) || isempty(EEG.data)
    errordlg('Empty EEG dataset.');
    return
end
if isempty(EEG.chanlocs(1).labels)
    errordlg('No channel labels.');
    return
end

% Get parameters from GUI or command line
if nargin == 1
    [params, abort] = getparams_gui(EEG);       % GUI
    if abort
        disp('Aborted.'); return
    end
elseif nargin > 1
    params = getparams_cmd(varargin{:});    % Command line
end

% run routine checks
[EEG, params, err] = run_checks(EEG,params);
if err, return; end  % stop program if there was an error (because errdlg
% does not stop program)

% Keep only EEG and ECG channels in separate structures
% EEG = pop_select(EEG, 'chantype',{'ECG','EEG'});
CARDIO = pop_select(EEG,'channel',params.heart_channels); % export ECG data in separate structure
% if isempty(params.heart_channels)  % for computing EEG features without heart data
%     CARDIO.data = [];
% end
EEG = pop_select(EEG,'nochannel',params.heart_channels);


% Basic check for other auxiliary channels that could cause issues
% (very simplistic and limited by dataset's electrode labels)
idx = contains(lower({EEG.chanlocs.labels}), {'ecg' 'ekg' 'ppg' 'aux' 'gsr' 'eda' 'eog' 'emg'});
if any(idx)
    auxChan = strjoin({EEG.chanlocs(idx).labels});
    opts.Interpreter = 'tex';
    opts.Default = 'Yes';
    quest = sprintf("The following channels may not be EEG or ECG/PPG and may cause serious innacuracies: \n\n %s \n\n They are about to be removed, please confirm", auxChan);
    answer = questdlg(quest,'Potentially undesirable electrodes','Yes','No','Cancel',opts);
    % ans = questdlg(");
    if strcmp(answer,'Cancel')
        disp('Aborted.'); return
    elseif strcmp(answer,'Yes')
        EEG = pop_select(EEG,'nochannel',{EEG.chanlocs(idx).labels});
    end
end

%%%%% MODE 1 & 2: RR, SQI, and NN %%%%%
if ~strcmpi(params.heart_signal,'off')
    if strcmp(params.analysis, 'hep') || isfield(params,'hrv_features')

        % Resample CARDIO data to match EEG for HEP, when different
        if  strcmp(params.analysis, 'hep') && EEG.srate~=CARDIO.srate
            fprintf('Resampling cardiovascular data to match EEG sampling rate. \n')
            CARDIO = pop_resample(CARDIO,EEG.srate);
        end

        % Get RR and NN intervals from ECG/PPG signals
        % note: when several electrodes are provided, use the elec with best
        % quality for subsequent analysis. Structures are used to avoid issues
        % when outputs have different lengths across electrodes.
        signal = CARDIO.data;
        sqi = [];
        nElec = size(signal,1);
        for iElec = 1:nElec
            elec = sprintf('elec%g',iElec);
            fprintf('Detecting R peaks from cardiovascular time series %g (%s)... \n', iElec, CARDIO.chanlocs(iElec).labels)
            [RR.(elec), RR_t.(elec), Rpeaks.(elec), sig(iElec,:), sig_t(iElec,:), pol.(elec), HR(iElec,:)] = get_RR(signal(iElec,:)',params);

            % % Fix values if PPG had a different sampling rate than EEG (this
            % should now be avoided by resampling above)
            factor = EEG.srate / CARDIO.srate;
            if factor ~= 1
                errordlg("Your EEG and cardiovascular data must have the same sampling rate.")
                %     RR.(elec) = RR.(elec) .* factor;
                %     RR_t.(elec) = (Rpeaks.(elec) .* factor) / EEG.srate;
                %     % RR_t.(elec) = (Rpeaks(1:end-1) .* factor) / EEG.srate;
                %     CARDIO = pop_resample(CARDIO,EEG.srate);
                %     % sig = CARDIO.data(iElec,:);
                %     sig_t(iElec,:) = sig_t(iElec,:) .* factor;
                %     params.fs = CARDIO.srate;
            end

            % SQI
            SQIthresh = .9; % minimum SQI recommended by Vest et al. (2019)
            maxThresh = 20; % max portion of artifacts (20% default from Vest et al. 2019)
            if strcmpi(params.heart_signal, 'ecg')
                sqi(iElec,:) = get_sqi_ecg(Rpeaks.(elec), signal(iElec,:), params.fs);
            elseif strcmpi(params.heart_signal, 'ppg')
                [sqi(iElec,:),~,annot] = get_sqi_ppg(Rpeaks.(elec),signal(iElec,:)',params.fs);
            else
                errdlg("Heart channel should be either 'ECG' or 'PPG' ")
            end
            SQI_mu(iElec,:) = round(mean(sqi(iElec,:), 'omitnan'),2);
            SQI_badRatio(iElec,:) = round(sum(sqi(iElec,:) < SQIthresh) / length(sqi(iElec,:))*100,1);
            % if SQI_mu(iElec,:) < .9
            %     warning("Mean signal quality index (SQI): %g. Minimum recommended SQI = .9 before correction of RR artifacts. See Vest et al. (2017) for more detail. \n", SQI_mu)
            % else
            % end
            % if SQI_badRatio(iElec,:) > 20
            %     warning("%g%% of the signal quality index (SQI) is below the minimum threshold (SQI = .9) before correction of RR artifacts. Maximum recommended is %g%%. See Vest et al. (2017) for more detail. \n", SQI_badRatio(iElec,:), maxThresh)
            %     warndlg(sprintf("%g%% of the signal quality index (SQI) is below the minimum threshold (SQI = .9) before correction of RR artifacts. Maximum portion recommended is %g%%. See Vest et al. (2017) for more detail. \n", SQI_badRatio(iElec,:), maxThresh),'Signal quality warning 1')
            % else
            %     fprintf("%g%% of the signal quality index (SQI) is below the minimum threshold (SQI = .9) before correction of RR artifacts. Maximum recommended is %g%%. See Vest et al. (2017) for more detail. \n", SQI_badRatio(iElec,:), maxThresh)
            % end

            % Correct RR artifacts (e.g., arrhytmia, ectopy, noise) to obtain the NN series
            disp("Correcting RR artifacts...")
            warning off
            [NN.(elec), NN_t.(elec), flagged.(elec)] = clean_rr(RR_t.(elec), RR.(elec), params);
            flaggedRatio.(elec) = sum(flagged.(elec)) / length(flagged.(elec)) *100;
            warning on

        end

        % Keep only ECG data of electrode with the lowest number of RR
        % artifacts
        [~,best_elec] = min(struct2array(flaggedRatio));
        elec = sprintf('elec%g',best_elec);
        flaggedRatio = flaggedRatio.(elec);
        flagged = flagged.(elec);
        if sum(flagged) > 0
            fprintf('%g/%g (%g%%) of heart beats were flagged as artifacts and interpolated (or removed if you chose to remove them). \n', sum(flagged),length(flagged),round(flaggedRatio,2));
        end
        if flaggedRatio > maxThresh % more than 20% of RR series is bad file
            warning("%g%% of the RR series on your best electrode was flagged as artifact. Maximum recommendation is 20%%. You may want to check for abnormal sections (e.g. electrode disconnections for long periods of time) in your cardiovascular signal and try BrainBeats again. ", round(flaggedRatio,2));
            warndlg(sprintf("%g%% of the RR series on your best electrode was flagged as artifact. Maximum recommendation is 20%%. You may want to check for abnormal sections (e.g. electrode disconnections for long periods of time) in your cardiovascular signal and try BrainBeats again.", round(flaggedRatio,2)),'Signal quality warning 2');
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
        pol = pol.(elec); % ECG signal polarity
        SQI_mu = SQI_mu(best_elec,:);
        SQI_badRatio = SQI_badRatio(best_elec,:);

        % Print average SQI
        fprintf("Overall signal quality index (SQI): %g. Note: SQI > 0.9 is considered good.\n", SQI_mu)

        % Plot filtered ECG and RR series of best electrode and interpolated
        % RR artifacts (if any)
        if params.vis_cleaning
            plot_NN(sig_t, sig, RR_t, RR, Rpeaks, NN_t, NN, flagged, params.heart_signal)
            pause(0.1)  % to avoid waiting for EEG preprocessing to appear
        end

        % Preprocessing outputs
        % EEG.brainbeats.preprocessing.RR_times = RR_t;
        % EEG.brainbeats.preprocessing.RR = RR;
        EEG.brainbeats.preprocessings.interpolated_heartbeats = sum(flagged);
        if exist('sqi','var')
            EEG.brainbeats.preprocessings.heart_SQI_mean = SQI_mu;
            EEG.brainbeats.preprocessings.heart_SQI_badportion = SQI_badRatio;
        end
        EEG.brainbeats.preprocessings.NN = NN;
        EEG.brainbeats.preprocessings.NN_times = NN_t;

        % % Take filtered cardio signal
        % if strcmp(params.heart_signal,'ecg') && pol<0
        %     CARDIO.data = -sig; % reverse to positive polarity
        % else
        %     CARDIO.data = sig;
        % end

        %%%%% MODE 2: HRV features %%%%%
        if ~isempty(params.hrv_features) && params.hrv_features ~= 0

            % File length (we take the whole series for now to allow ULF and VLF as much as possible)
            file_length = floor(EEG.xmax)-1;
            if file_length < 300
                warning('File length is less than 5 minutes! The minimum recommended is 300 s for estimating reliable HRV metrics.')
                warndlg('File length is less than 5 minutes! The minimum recommended is 300 s for estimating reliable HRV metrics.')
            end
            params.file_length = file_length;

            % Extract HRV measures
            [features_hrv, params] = get_hrv_features(NN, NN_t, params);

            % Final output with everything
            Features.HRV = features_hrv;
            Features.HRV.time.heart_rate = round(mean(HR,'omitnan'),1);

            % Exit BrainBeats if user only wants to work with cardiovascular data
            if strcmp(params.analysis,'features') && ~params.eeg_features
                EEG.brainbeats.features = Features;
                disp('Done processing cardiovascular signals'); gong
                % return
            end
        end
    end
end

% Basic preprocessing of EEG signals (filter, reference, bad channels)
if params.clean_eeg
    params.clean_eeg_step = 0;
    params.orichanlocs = EEG.chanlocs;
    [EEG, params] = clean_eeg(EEG, params);

    % Preprocessing outputs
    EEG.brainbeats.preprocessings.removed_eeg_channels = params.removed_eeg_channels;
end

%%%%% Heartbeat-evoked potentials (HEP) %%%%%
if strcmp(params.analysis,'hep')
    EEG = run_HEP(EEG, CARDIO, params, Rpeaks); % CARDIO input is for adding cardio channel back in final output
end

%%%%% EEG features %%%%%
if isfield(params,'eeg_features') && ~isempty(params.eeg_features) && params.eeg_features~=0

    % Remove EEG artifacts with ASR and bad components with ICLabel
    if params.clean_eeg
        [EEG, params] = clean_eeg(EEG, params);

        % Preprocessing outputs
        EEG.brainbeats.preprocessings.removed_eeg_segments = params.removed_eeg_segments;
        EEG.brainbeats.preprocessings.removed_eeg_components = params.removed_eeg_components;
    end

    % Extract EEG features and store in EEGLAB structure
    params.chanlocs = EEG.chanlocs;
    [features_eeg, params] = get_eeg_features(EEG.data, params);

    % Final output with everything
    Features.EEG = features_eeg;
end

% Store parameters in EEG structure for reporting
EEG.brainbeats.parameters = params;

% Store, save, plot features
if strcmp(params.analysis,'features')
    EEG.brainbeats.features = Features;
    disp("Features are stored in the EEG.features field")

    % Save in same repo as original file
    if params.save
        outputPath = fullfile(EEG.filepath, sprintf('%s_features.mat', EEG.filename(1:end-4)));
        fprintf("Exporting all features in a .mat file at this location: %s \n", outputPath);
        save(outputPath,'Features');
    end

    % Plot features
    if params.vis_outputs
        plot_features(Features,params)
    end
end

%%%%% MODE 3: remove heart components from EEG signals %%%%%
if strcmp(params.analysis,'rm_heart')

    % Preprocess EEG
    if params.clean_eeg

        % ref, filt, bad channels
        params.clean_eeg_step = 0;
        [EEG, params] = clean_eeg(EEG, params);

        % EEG artifacts with ASR and ICA
        params.clean_eeg_step = 1;
        [EEG, params] = clean_eeg(EEG, params);

        % Filter ECG
        if ~isfield(params,'highpass_ecg')
            params.highpass_ecg = 1;
        end
        if ~isfield(params,'lowpass_ecg')
            params.lowpass_ecg = 20;
        end
        CARDIO = pop_eegfiltnew(CARDIO,'locutoff',params.highpass_ecg,'hicutoff',params.lowpass_ecg);

        % Remove same segments from cardio signal
        CARDIO = pop_select(CARDIO,'nopoint', params.removed_eeg_segments);

        % Preprocessing outputs
        EEG.brainbeats.preprocessings.removed_eeg_segments = params.removed_eeg_segments;
        EEG.brainbeats.preprocessings.removed_eeg_components = params.removed_eeg_components;
    end

    % Add CARDIO channel back with EEG
    EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
    EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
    for iChan = 1:CARDIO.nbchan
        EEG.chanlocs(end+1).labels = params.heart_channels{iChan};
    end
    EEG = eeg_checkset(EEG);

    % if params.vis_cleaning
    %     pop_eegplot(EEG,1,1,1);
    % end

    EEG = remove_heartcomp(EEG, params);
end


% Command history (for 'eegh')
if contains(params.analysis,'hep')
    com = char(sprintf("EEG = brainbeats_process(EEG,'analysis','%s','heart_signal','%s','heart_channels',%s,'clean_eeg',%g,'parpool',%g,'gpu',%g,'vis_cleaning',%g,'vis_outputs',%g,'save',%g);",...
        params.analysis,params.heart_signal,['{' sprintf(' ''%s'' ', params.heart_channels{:}) '}'],params.parpool,params.gpu,params.vis_cleaning,params.vis_outputs,params.save));
elseif contains(params.analysis,'features')
    if params.hrv_features
        opt = {'time' 'frequency' 'nonlinear'};
        idx = logical([params.hrv_time params.hrv_frequency params.hrv_nonlinear]);
        hrv_features = ['{' sprintf(' ''%s'' ', opt{idx}) '}'];
    else
        hrv_features = ' ''off'' ';
    end
    if params.eeg_features
        opt = {'time' 'frequency' 'nonlinear'};
        idx = logical([params.eeg_time params.eeg_frequency params.eeg_nonlinear]);
        eeg_features = ['{' sprintf(' ''%s'' ', opt{idx}) '}'];
    else
        eeg_features = ' ''off'' ';
    end
    com = char(sprintf("EEG = brainbeats_process(EEG,'analysis','%s','heart_signal','%s','heart_channels',%s,'clean_eeg',%g,'hrv_features',%s,'eeg_features',%s,'parpool',%g,'gpu',%g,'vis_cleaning',%g,'vis_outputs',%g,'save',%g);",...
        params.analysis,params.heart_signal,['{' sprintf(' ''%s'' ', params.heart_channels{:}) '}'],params.clean_eeg,hrv_features,eeg_features,params.parpool,params.gpu,params.vis_cleaning,params.vis_outputs,params.save));
elseif contains(params.analysis,'rm_heart')
    com = char(sprintf("EEG = brainbeats_process(EEG,'analysis','%s','heart_signal','%s','heart_channels',%s,'clean_eeg',%g,'vis_cleaning',%g,'vis_outputs',%g,'save',%g);",...
        params.analysis,params.heart_signal,['{' sprintf(' ''%s'' ', params.heart_channels{:}) '}'],params.clean_eeg,params.vis_cleaning,params.vis_outputs,params.save));
end

% Final message with ref to cite
fprintf('\n')
fprintf("Done! Thank you for using the BrainBeats toolbox! Please cite: \n");
fprintf("Cannard, C., Wahbeh, H., Delorme, A. BrainBeats as an Open-Source EEGLAB Plugin to Jointly Analyze EEG and Cardiovascular Signals. J. Vis. Exp. (2024). \n")
fprintf("https://review.jove.com/t/65829/brainbeats-as-an-open-source-eeglab-plugin-to-jointly-analyze-eeg?status=a67835k \n")

if params.gong
    gong
end
