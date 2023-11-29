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
%  'vis_cleaning'   - [0|1] set vizualization of preprocessing plots to 
%                     on (1) or off (0). Default is on.
%  'vis_outputs'    - [0|1] set vizualization of outputs to on (1) or 
%                     off (0). Default is on.
%  'save'           - [0|1] save results into a MATLAB file (1) or not (0). Default is on.
%
% Outputs:
%   EEG      - modified EEGLAB dataset. A EEG.features field is created
%              containing all the features computed for the dataset.
%
% Copyright (C) - BrainBEats, Cedric Cannard, 2023

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

% Get parameters from GUI or command line
if nargin == 1
    [params, abort] = getparams_gui(EEG);                % GUI
    if abort
        disp('Aborted.'); return 
    end
elseif nargin > 1
    params = getparams_command(varargin{:});    % Command line
end

% run routine checks
[EEG, params, err] = run_checks(EEG,params);
if err, return; end  % stop program if there was an error (because errdlg 
% does not stop program)

%%%%% MODE 1 & 2: RR, SQI, and NN %%%%%
if contains(params.analysis, {'features' 'hep'})
    
    % Keep only EEG and ECG channels in separate structures
    % EEG = pop_select(EEG, 'chantype',{'ECG','EEG'});  
    CARDIO = pop_select(EEG,'channel',params.heart_channels); % export ECG data in separate structure
    EEG = pop_select(EEG,'nochannel',params.heart_channels); 

    % basic check for missed auxiliary channels that could cause issues 
    % (very limited by dataset's electrode labels)
    idx = contains({EEG.chanlocs.labels}, {'ecg' 'ECG' 'aux' 'AUX'});
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
        [NN.(elec), NN_t.(elec), flagged.(elec)] = clean_rr(RR_t.(elec), RR.(elec), params);
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
    %     errordlg("Signal quality is too low: aborting! You could inspect the data in EEGLAB > Plot > Channel data (Scroll) and try to remove large artifacts first.");
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
    if params.vis_cleaning
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

    %%%%% MODE 1: Heartbeat-evoked potentials (HEP) %%%%%
    if strcmp(params.analysis,'hep')
        EEG = run_HEP(EEG, params, Rpeaks);
    end 
    
    %%%%% MODE 2: HRV features %%%%%
    if strcmp(params.analysis,'features') && params.hrv

            % File length (we take the whole series for now to allow ULF and VLF as much as possible)
            file_length = floor(EEG.xmax)-1;
            if file_length < 300
                warning('File length is shorter than 5 minutes! The minimum recommended is 300 s for estimating reliable HRV metrics.')
            end
            params.file_length = file_length;

            % Extract HRV measures
            features_hrv = get_hrv_features(NN, NN_t, params);

            % Final output with everything
            Features.HRV = features_hrv;
    end
    
    %%%%% MODE 2: EEG features %%%%%
    if strcmp(params.analysis,'features') && params.eeg

        % Clean EEG artifacts with ASR
        if params.clean_eeg
            [EEG, params] = clean_eeg(EEG, params);
        end

        % Extract EEG features and store in EEGLAB structure
        if params.eeg
            params.chanlocs = EEG.chanlocs;
            features_eeg = get_eeg_features(EEG.data, params);
        end

        % Final output with everything
        Features.EEG = features_eeg;

    end

    % Store in EEGLAB structure (independent of saving choice in case user 
    % wants to save .set file on their own in a script)
    if strcmp(params.analysis,'features')
        EEG.features = Features;
        disp("Features are stored in the EEG.features field")
    end

    % Save and plot features  
    if strcmp(params.analysis,'features')

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

end

%%%%% MODE 3: remove heart components from EEG signals %%%%%
if strcmp(params.analysis,'rm_heart')
    EEG = remove_heartcomp(EEG, params);
end

% Shut down parallel pool
% if params.parpool
%     delete(gcp('nocreate'));
% end


% eegh
if contains(params.analysis,'hep')
    com = char(sprintf("EEG = brainbeats_process(EEG,'analysis','%s','heart_signal','%s','heart_channels','%s','rr_correct','%s','clean_eeg',%g,'parpool',%g,'gpu',%g,'vis_cleaning',%g,'vis_outputs',%g,'save',%g);",...
        params.analysis,params.heart_signal,['{' sprintf(' ''%s'' ', params.heart_channels{:}) '}'],params.rr_correct,params.clean_eeg,params.parpool,params.gpu,params.vis_cleaning,params.vis_outputs,params.save));
elseif contains(params.analysis,'features') 
    com = char(sprintf("[~, Features] = brainbeats_process(EEG,'analysis','%s','heart_signal','%s','heart_channels','%s','rr_correct','%s','clean_eeg',%g,'parpool',%g,'gpu',%g,'vis_cleaning',%g,'vis_outputs',%g,'save',%g);",...
        params.analysis,params.heart_signal,['{' sprintf(' ''%s'' ', params.heart_channels{:}) '}'],params.rr_correct,params.clean_eeg,params.parpool,params.gpu,params.vis_cleaning,params.vis_outputs,params.save));
elseif contains(params.analysis,'rm_heart') 
    com = char(sprintf("EEG = brainbeats_process(EEG,'analysis','%s','heart_signal','%s','heart_channels','%s','clean_eeg',%g,'vis_cleaning',%g,'vis_outputs',%g,'save',%g);",...
        params.analysis,params.heart_signal,['{' sprintf(' ''%s'' ', params.heart_channels{:}) '}'],params.clean_eeg,params.vis_cleaning,params.vis_outputs,params.save));
end
com = char(com);

% References
% if params.hrv
%     fprintf("For cardivosacular signal processing and HRV metrics, please cite: \n ");
%     fprintf("   - Vest et al. (2018). An open source benchmarked toolbox for cardiovascular waveform and interval analysis. Physiol Meas. \n");
%     fprintf("   - Shaffer & Ginsberg (2017). An overview of heart rate variability metrics and norms. Frontiers in public health. \n");
% end
% if strcmp(params.analysis,'features')
%     if params.eeg && params.eeg_frequency
%         fprintf("For the IAF feature, please cite: \n %s \n", "  - Corcoran et al. (2018). Toward a reliable, automated method of individual alpha frequency (IAF) quantification. Psychophysiology. ")
%         fprintf("For the alpha asymmetry feature, please cite: \n %s \n", "  - Smith et al. (2017). Assessing and conceptualizing frontal EEG asymmetry: An updated primer on recording, processing, analyzing, and interpreting frontal alpha asymmetry. International Journal of Psychophysiology.")
%         fprintf("For the EEG coherence measures, please cite: \n");
%         fprintf("   - Cao et al. (2016). Resting-state EEG power and coherence vary between migraine phases. The journal of headache and pain. \n");
%         fprintf("   - Pascual-Marqui et al. (2014). Assessing direct paths of intracortical causal information flow of oscillatory activity with the isolated effective coherence (iCoh). Frontiers in human neuroscience. \n");
%     end
%     if params.hrv_nonlinear || params.eeg_nonlinear
%         fprintf("For the entropy measures, please cite: \n");
%         fprintf("   - Azami & Escudero (2016). Refined Multiscale Fuzzy Entropy based on standard deviation for biomedical signal analysis. Medical & Biological Engineering & Computing. \n");
%         fprintf("   - Costa, Goldberger, Peng (2002). Multiscale entropy analysis of complex physiologic time series. Phys Rev Lett. \n ")
%         fprintf("   - Kosciessa et al. (2020). Standard multiscale entropy reflects neural dynamics at mismatched temporal scales: What's signal irregularity got to do with it? Plos Comput Biol. \n ")
%     end
% end
% if strcmp(params.analysis,'hep')
%     fprintf("For HEP methods, please cite: \n");
%     fprintf("   - Candia-Rivera et al. (2021). The role of electroencephalography electrical reference in the assessment of functional brain–heart interplay: From methodology to user guidelines. Neuroscience Methods. \n");
%     fprintf("   - Park & Blanke (2019). Heartbeat-evoked cortical responses: Underlying mechanisms, functional roles, and methodological considerations. NeuroImage. \n");
% end
fprintf('\n\n')
fprintf("Done! Thank you for using the BrainBeats toolbox! \n");
fprintf("Please cite: Cannard, Wahbeh, & Delorme (2023). BrainBeats: an open-source EEGLAB plugin to jointly analyze EEG and cardiovascular (ECG/PPG) signals. bioRxiv, 2023-06 \n")

 gong
