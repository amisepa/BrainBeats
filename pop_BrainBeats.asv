%% Main script
%
% Potential names: BrainBeats, CardioNeuroSync (CNS), NeuroPulse, CardioCortex
%
% INPUTS:
%   'rr_correct' - method to interpolate RR artifacts (e.g. 'linear','cubic')
%
% Cedric Cannad, 2023

function [EEG, Features, com] = pop_BrainBeats(EEG, varargin)

pop_editoptions('option_single', 0); % ensure double precision
Features = [];
com = '';

% Add path to subfolders
mainpath = fileparts(which('pop_BrainBeats.m'));
addpath(fullfile(mainpath, 'functions'));
addpath(fullfile(mainpath, 'functions', 'restingIAF'));

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

%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%
if nargin == 1
    params = getparams_gui(EEG);                % GUI
elseif nargin > 1
    params = getparams_command(varargin{:});    % Command line
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if GUI was aborted
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
end

% Check for channel locations for visualization
if params.vis
    if ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(1).X)
        error("Electrode location coordinates must be loaded for visualizing outputs.");
    end
end

% GPU computing
if ~isfield(params,'gpu') % not available from GUI yet
    params.gpu = false;
end

%%%%%%%%%%%%% PREPROCESS EEG DATA %%%%%%%%%%%%%

EEG.data = double(EEG.data);  % ensure double precision
params.fs = EEG.srate;
if ~iscell(params.heart_channels)
    params.heart_channels = {params.heart_channels};
end
ECG = pop_select(EEG,'channel',params.heart_channels); % export ECG data in separate structure

% Filter, re-reference, remove bad channels
if params.clean_eeg
    params.clean_eeg_step = 0;
    [EEG, params] = clean_eeg(EEG, params);
end


%%%%% MODE 1: remove heart components from EEG signals with IClabel %%%%%

if strcmp(params.analysis,'rm_heart')
    remove_heartcomp(EEG, params);
end

%%%%% MODE 2 & 3: RR, SQI, and NN %%%%%
if contains(params.analysis, {'features' 'hep'})

    % Get RR series and signal quality index (SQI)
    if strcmp(params.heart_signal, 'ecg')

        % idx = contains({EEG.chanlocs.labels}, params.heart_channels);
        % ecg = EEG.data(idx,:);
        % ECG = pop_resample(ECG,125);  % For get_rwave2
        ecg = ECG.data;
        nElec = size(ecg,1);
        for iElec = 1:nElec
            elec = sprintf('elec%g',iElec);
            fprintf('Detecting R peaks from ECG time series: electrode %g...\n', iElec)
            [RR.(elec), RR_t.(elec), Rpeaks.(elec), sig_filt(iElec,:), sig_t(iElec,:), HR] = get_RR(ecg(iElec,:)', params);

            % SQI
            SQIthresh = .9;
            [sqi(iElec,:), sqi_times(iElec,:)] = get_sqi(Rpeaks.(elec), ecg(iElec,:), params.fs);
            SQI(iElec,:) = sum(sqi(iElec,:) < SQIthresh) / length(sqi(iElec,:));  % minimum SQI recommended by Vest et al. (2019)
        end

        % Keep only ECG data of electrode with the best SQI
        [~, best_elec] = min(SQI);
        % sqi = sqi(best_elec,:);
        sqi = [sqi_times(best_elec,:); sqi(best_elec,:)];

        SQI = SQI(best_elec);
        sig_t = sig_t(best_elec,:);
        sig_filt = sig_filt(best_elec,:);
        SQIthresh2 = .2;   % 20% of file can contain SQI<.9
        if SQI > SQIthresh2 % more than 20% of RR series is bad
            warning(['%g%% of the RR series on your best ECG electrode has a signal quality index (SQI) below minimum recommendations (max 20%% below SQI = .9; see Vest et al., 2019)! \n' ...
                'You may inspect and remove them manually in EEGLAB > Plot > Channel data (Scroll).'], SQI)
        else
            fprintf( "Keeping only the heart electrode with the best signal quality index (SQI): %g%% of the RR series is outside of the recommended threshold. \n", SQI )
        end
        elec = sprintf('elec%g',best_elec);
        RR = RR.(elec);
        RR_t = RR_t.(elec);
        RR_t(1) = [];       % always ignore 1st hearbeat
        Rpeaks = Rpeaks.(elec);
        Rpeaks(1) = [];     % always ignore 1st hearbeat

    elseif strcmp(params.heart_signal,'ppg')
        error("Work in progress, sorry!");
        % [rr,t_rr,sqi] = Analyze_ABP_PPG_Waveforms(InputSig,{'PPG'},HRVparams,[],subID);

    else
        error("Unknown heart signal. Should be 'ecg' or 'ppg' ");
    end

    % Correct RR artifacts (e.g., arrhytmia, ectopy, noise) to obtain the NN series
    vis = false;    % to visualize artifacts that are inteprolated
    [NN, NN_times,flagged] = clean_rr(RR_t, RR, sqi, params, vis);
    if sum(flagged) > 0
        if contains(params.rr_correct,'remove')
            fprintf('%g heart beats were flagged as artifacts and removed. \n', sum(flagged));
        else
            fprintf('%g heart beats were flagged as artifacts and interpolated. \n', sum(flagged));
        end
    end

    % Outputs
    Features.HRV.ECG_filtered = sig_filt;
    Features.HRV.ECG_times = sig_t;
    Features.HRV.SQI = SQI;
    Features.HRV.RR = RR;
    Features.HRV.RR_times = RR_t;
    Features.HRV.HR = HR;
    Features.HRV.NN = NN;
    Features.HRV.NN_times = NN_times;
    Features.HRV.flagged_heartbeats = flagged;

    % Plot filtered ECG and RR series of best electrode and interpolated
    % RR artifacts (if any)
    if params.vis
        figure('color','w');

        subplot(2,1,1)
        % scrollplot({sig_t,sig_filt,'color','#0072BD'},{RR_t,sig_filt(Rpeaks),'.','MarkerSize',10,'color','#D95319'}, {'X'},{''},.2);
        plot(sig_t, sig_filt,'color','#0072BD'); hold on;
        plot(RR_t, sig_filt(Rpeaks),'.','MarkerSize',10,'color','#D95319');
        % title(sprintf('Filtered ECG signal + R peaks (portion of artifacts: %1.2f%%)',SQI)); ylabel('mV'); %set(gca,'XTick',[]);
        title('Filtered ECG signal + R peaks'); ylabel('mV'); axis tight %set(gca,'XTick',[]);

        subplot(2,1,2)
        if sum(flagged) == 0
            plot(RR_t,RR,'-','color','#0072BD','linewidth',1);
        else
            plot(RR_t,RR,'-','color','#A2142F','linewidth',1);
            hold on; plot(NN_times, NN,'-','color',"#0072BD", 'LineWidth', 1);
            legend('RR artifacts','NN intervals')
        end
        title('RR intervals'); ylabel('RR intervals (s)'); xlabel('Time (s)'); axis tight
        set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold'); box on
    end

    %%%%% MODE 2: Heartbeat-evoked potentials (HEP) %%%%%
    if strcmp(params.analysis,'hep')
        EEG = run_HEP(EEG, params, Rpeaks);
    end

    %%%%% MODE 3: HRV features %%%%%
    if strcmp(params.analysis,'features') && params.hrv

        if SQI <= .2 % <20% of RR interval is artifacts

            % defaults
            params.hrv_norm = true;  % default
            params.hrv_spec = 'Lomb-Scargle periodogram';  % 'Lomb-Scargle periodogram' (default), 'pwelch', 'fft', 'burg'
            params.hrv_overlap =  0.25; % 25%

            % File length (we take the whole series for now to allow ULF and VLF as much as possible)
            file_length = floor(EEG.xmax)-1;
            if file_length < 300
                warning('File length is shorter than 5 minutes! The minimum recommended is 300 s for estimating reliable HRV metrics.')
            end
            params.file_length = file_length;

            % Extract HRV measures
            HRV = get_hrv_features(NN, NN_times, params);

            % Final output with everything
            Features.HRV = HRV;
            
        else
            error('Signal quality of the RR series is too low. HRV features should not be computed.')
        end
    end

    %%%%% MODE 3: EEG features %%%%%
    if strcmp(params.analysis,'features') && params.eeg

        % Clean EEG artifacts with ASR
        if params.clean_eeg
            [EEG, params] = clean_eeg(EEG, params);
        end

        % Extract EEG features
        if params.eeg
            tic
            eeg_features = get_eeg_features(EEG.data, EEG.chanlocs, params);
            toc
        end

        % Final output with everything
        Features.EEG = eeg_features;

    end


    %%%%%% PLOT EEG AND HRV FEATURES %%%%%%%
    
    % Visualize HRV outputs
    if strcmp(params.analysis,'features') && params.vis
        plot_features(Features)
    end

end

save(fullfile(outDir, 'features.mat'),'Features'); %FIXME: ASK USER FOR OUTPUT DIR

disp('Done!'); gong


