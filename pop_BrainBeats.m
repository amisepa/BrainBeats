%% Main script
%
% Potential names: BrainBeats, CardioNeuroSync (CNS), NeuroPulse, CardioCortex
%
% Cedric Cannad, 2023

function [outputs, com] = pop_BrainBeats(EEG, varargin)

outputs = [];
com = '';

% Add path to subfolders
mainpath = fileparts(which('pop_BrainBeats.m'));
addpath(fullfile(mainpath, 'functions'));

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

%%%%%%%%%%%%% PREPROCESS DATA %%%%%%%%%%%%%

EEG.data = double(EEG.data);  % ensure double precision
ECG = pop_select(EEG,'channel',params.heart_channels); % export ECG data in separate structure
if params.clean_eeg
    EEG = pop_eegfiltnew(EEG,'locutoff',1,'hicutoff',45,'filtorder',846);
    EEG = pop_select(EEG,'nochannel',params.heart_channels); % FIXME: remove all non-EEG channels instead

    % Reference to infinity
    if ~isfield(EEG,'ref') || isempty(EEG.ref) || strcmp(EEG.ref,'')
        EEG = reref_inf(EEG); % my function
    end

    % Remove bad channels
    oriEEG = EEG;
    EEG = pop_clean_rawdata(EEG,'FlatlineCriterion',10,'ChannelCriterion',.85, ...
        'LineNoiseCriterion',5,'Highpass','off', 'BurstCriterion','off', ...
        'WindowCriterion','off','BurstRejection','off','Distance','off');

    % Identify periods with large artifacts using ASR
    % cutoff = 60;
    % useriemannian = false;
    % m = memory;
    % maxmem = round(.85*(m.MemAvailableAllArrays/1000000),1);  % use 85% of available memory (in MB)
    % cleanEEG = clean_asr(EEG,cutoff,[],[],[],[],[],[],false,useriemannian,maxmem);
    % mask = sum(abs(EEG.data-cleanEEG.data),1) > 1e-10;
    % EEG.etc.clean_sample_mask = ~mask;
    % badData = reshape(find(diff([false mask false])),2,[])';
    % badData(:,2) = badData(:,2)-1;
    % if ~isempty(badData) % ignore very small artifacts (<5 samples)
    %     smallIntervals = diff(badData')' < 5;
    %     badData(smallIntervals,:) = [];
    % end

    % % Remove bad segments
    % EEG = pop_select(EEG,'nopoint',badData);
    % if strcmp(params.analysis,'hep')
    %     ECG = pop_select(ECG,'nopoint',badData);
    % end
    % fprintf('%g %% of data were considered to be large artifacts and removed. \n', (1-EEG.xmax/oriEEG.xmax)*100)
    %
    % % Visualize what was removed
    if params.vis
        vis_artifacts(EEG,oriEEG,'ChannelSubset',1:EEG.nbchan-length(params.heart_channels));
    end

    % Interpolate bad channels
    EEG = pop_interp(EEG, oriEEG.chanlocs, 'spherical'); % interpolate
    EEG.etc.clean_channel_mask(1:EEG.nbchan) = true;

    % % Add ECG channels back
    % EEG.data(end+1:end+ECG.nbchan,:) = ECG.data;
    % EEG.nbchan = EEG.nbchan + ECG.nbchan;
    % for iChan = 1:ECG.nbchan
    %     EEG.chanlocs(end+1).labels = params.heart_channels{iChan};
    % end
    % EEG = eeg_checkset(EEG);

end

params.fs = EEG.srate;

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
    outputs.HRV.ECG_filtered = sig_filt;
    outputs.HRV.ECG_times = sig_t;
    outputs.HRV.SQI = SQI;
    outputs.HRV.RR = RR;
    outputs.HRV.RR_times = RR_t;
    outputs.HRV.HR = HR;
    outputs.HRV.NN = NN;
    outputs.HRV.NN_times = NN_times;
    outputs.HRV.flagged_heartbeats = flagged;

    % Plot filtered ECG and RR series of best electrode and interpolated
    % RR artifacts (if any)
    if params.vis
        figure('color','w');

        subplot(2,1,1)
        scrollplot({sig_t,sig_filt,'color','#0072BD'},{RR_t,sig_filt(Rpeaks),'.','MarkerSize',10,'color','#D95319'}, {'X'},{''},.2);
        % plot(sig_t, sig_filt,'color','#0072BD'); hold on;
        % plot(RR_t, sig_filt(Rpeaks),'.','MarkerSize',10,'color','#D95319');
        % title(sprintf('Filtered ECG signal + R peaks (portion of artifacts: %1.2f%%)',SQI)); ylabel('mV'); %set(gca,'XTick',[]);
        title('Filtered ECG signal + R peaks'); ylabel('mV'); %set(gca,'XTick',[]);

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
        run_HEP(EEG, params, Rpeaks)
    end

    %%%%% MODE 3: HRV features %%%%%
    if strcmp(params.analysis,'features') && params.hrv
        if SQI >= .2
            params.hrv_norm = true;  % default for now
            file_length = floor(EEG.xmax)-1;
            if file_length < 300
                warning('File length is shorter than 5 minutes! The minimum recommended is 300 s for estimating reliable HRV metrics.')
            end
            params.file_length = file_length;

            outputs = get_hrv_features(NN, NN_times, params);

        else
            error('Signal quality of the RR series is too low. HRV features should not be computed.')
        end
    end


    %%%%% MODE 3: EEG features %%%%%
    if strcmp(params.analysis,'features') && params.eeg

    end

end

disp('Done.')


