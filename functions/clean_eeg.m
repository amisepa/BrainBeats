function [EEG, params] = clean_eeg(EEG, params)

% Filter, re-reference, and remove bad channels
if params.clean_eeg_step == 0
    
    EEG = pop_eegfiltnew(EEG,'locutoff',1,'minphase',false);    % use causal minimum-phase filter for pre-event analysis
    EEG = pop_eegfiltnew(EEG,'hicutoff',45,'filtorder',846,'minphase',false); % use causal minimum-phase filter for pre-event analysis
    
    % Reference to average or infinity/REST
    % Candia-Rivera, Catrambone, & Valenza (2021). The role of EEG reference 
    % in the assessment of functional brain–heart interplay: From 
    % methodology to user guidelines. Journal of Neuroscience Methods.
    if isempty(EEG.ref) || strcmp(EEG.ref,'') || strcmp(EEG.ref,'unknown')
        if EEG.nbchan >= 30
            EEG = reref_inf(EEG); % my function
        else
            warning('Cannot reference these EEG data to infinity as low-density montages. Low-density montages require alternative referencing (e.g., linked mastoids).')
        end
    end
    
    % Remove bad channels
    oriEEG = EEG;
    EEG = pop_clean_rawdata(EEG,'FlatlineCriterion',5,'ChannelCriterion',.85, ...
        'LineNoiseCriterion',5,'Highpass','off', 'BurstCriterion','off', ...
        'WindowCriterion','off','BurstRejection','off','Distance','off');    
    % disp('Detecting flat line...')
    % EEG = clean_flatlines(EEG,5,20); 
    % try 
    %     EEG = clean_channels(EEG, .8, 4, 5, .6, 100); 
    % catch
    %     disp('Your dataset appears to lack correct channel locations; using a location-free channel cleaning method.');
    %     EEG = clean_channels_nolocs(EEG,.45,.2,10,.33,true); 
    % end
    badChan = ~contains({oriEEG.chanlocs.labels}, {EEG.chanlocs.labels});    
    fprintf(1, 'Bad channels were: ');
    fprintf(1, '%s ', oriEEG.chanlocs(badChan).labels )
    fprintf(1, '\n')
    
    % Visualize removed channels
    if params.vis
        % EEG.etc.clean_channel_mask(1:EEG.nbchan) = true;
        % EEG.etc.clean_channel_mask(badChan) = false;
        vis_artifacts(EEG,oriEEG);
        % vis_artifacts(EEG,oriEEG,'ChannelSubset',1:EEG.nbchan-length(params.heart_channels));
    end

    % Interpolate them
    EEG = pop_interp(EEG, oriEEG.chanlocs, 'spherical'); % interpolate
    EEG.etc.clean_channel_mask(1:EEG.nbchan) = true;

    % % Add ECG channels back
    % EEG.data(end+1:end+ECG.nbchan,:) = ECG.data;
    % EEG.nbchan = EEG.nbchan + ECG.nbchan;
    % for iChan = 1:ECG.nbchan
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

        % Remove bad trials
        warning('Removing bad trials...')
        b = design_fir(100,[2*[0 45 50]/EEG.srate 1],[1 1 0 0]);
        sigRMS = nan(1,size(EEG.data,3));
        snr = nan(1,size(EEG.data,3));
        for iEpoch = 1:size(EEG.data,3)
            sigRMS(:,iEpoch) = rms(rms(squeeze(EEG.data(:,:,iEpoch)),2));
            tmp = filtfilt_fast(b,1, squeeze(EEG.data(:,:,iEpoch))');
            snr(:,iEpoch) = rms(mad(squeeze(EEG.data(:,:,iEpoch)) - tmp'));
        end
        badRMS = isoutlier(sigRMS,'grubbs');
        badSNR = isoutlier(snr,'grubbs');
        badTrials = unique([find(badRMS) find(badSNR)]);
        % pop_eegplot(EEG,1,1,1);
        EEG = pop_rejepoch(EEG, badTrials, 0);

        % run RMS a 2nd time less aggressive in case some were missed
        sigRMS = nan(1,size(EEG.data,3));
        for iEpoch = 1:size(EEG.data,3)
            sigRMS(:,iEpoch) = rms(rms(squeeze(EEG.data(:,:,iEpoch)),2));
        end
        badRMS = isoutlier(sigRMS,'mean');
        EEG = pop_rejepoch(EEG, find(badRMS), 0);

    % Features
    elseif strcmp(params.analysis, 'features')

        % Identify artifacts using ASR
        oriEEG = EEG;
        cutoff = 40;
        useriemannian = false;
        m = memory;
        maxmem = round(.85*(m.MemAvailableAllArrays/1000000),1);  % use 85% of available memory (in MB)
        cleanEEG = clean_asr(EEG,cutoff,[],[],[],[],[],[],params.gpu,useriemannian,maxmem);
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
        %     ECG = pop_select(ECG,'nopoint',badData);
        % end
        fprintf('%g %% of data were considered to be artifacts and were removed. \n', (1-EEG.xmax/oriEEG.xmax)*100)
        
        if params.vis
            vis_artifacts(EEG,oriEEG); 
        end

    end
    
    % Run ICA at effective data rank (see Kim et al. 2023). Use Picard when
    % possible to increase speed.
    dataRank = sum(eig(cov(double(EEG.data(:,:)'))) > 1E-7);
    if exist('picard.m','file')
        EEG = pop_runica(EEG,'icatype','picard','maxiter',500,'mode','standard','pca',dataRank);
    else
        EEG = pop_runica(EEG,'icatype','runica','extended',1,'pca',dataRank);
    end

    % Classify and remove bad components with IClabel
    EEG = pop_iclabel(EEG,'default');
    EEG = pop_icflag(EEG,[NaN NaN; .95 1; .9 1; .95 1; NaN NaN; NaN NaN; NaN NaN]);
    badComp = find(EEG.reject.gcompreject);
    EEG = eeg_checkset(EEG);
    
    % Visualize tagged components
    if params.vis
        pop_selectcomps(EEG,1:20);         
        colormap("parula")
    end
    
    % Remove bad components
    if ~isempty(badComp)
        fprintf('Removing %g bad component(s). \n', length(badComp));
        EEG = pop_subcomp(EEG, badComp, 0);
    end
    
    % plot final cleaned data
    if params.vis
        if strcmp(params.analysis, 'hep')
            pop_eegplot(EEG,1,1,1);
        % elseif strcmp(params.analysis, 'features')
        %     eegplot(EEG.data,'winlength',15,'srate',EEG.srate,'events', ...
        %         EEG.event,'spacing',50);
        end
    end

    % update tracker
    params.clean_eeg_step = 2;

end 
