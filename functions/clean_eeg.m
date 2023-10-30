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
    EEG.etc.clean_channel_mask = true(1,EEG.nbchan);
    oriEEG = EEG;
    % EEG = pop_clean_rawdata(EEG,'FlatlineCriterion',5,'ChannelCriterion',.85, ...
    %     'LineNoiseCriterion',5,'Highpass','off', 'BurstCriterion','off', ...
    %     'WindowCriterion','off','BurstRejection','off','Distance','off');    
    disp('Detecting bad EEG channels...')
    EEG = clean_flatlines(EEG,5);   % remove channels that have flat lines
    win_length = 5;     % window length (deafult = 5 s)
    line_thresh = 5;    % line noise threshold (default = 5)
    maxBad = .15;        % max tolerated portion of channel to be bad before removal (default = .4)
    corrThresh = .85;   % correlation threshold to be considered bad (default = .85)
    nSamp = 200;        % number of ransac samples (default = 50; higher is more robust but longer)
    try 
        EEG = clean_channels(EEG,corrThresh,line_thresh,win_length,maxBad,nSamp); 
    catch
        warning('Your dataset appears to lack correct channel locations. Using the location-free parameters to remove bad EEG channels.');
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
    
    % Visualize removed channels
    if params.vis
        % EEG.etc.clean_channel_mask(1:EEG.nbchan) = true;
        % EEG.etc.clean_channel_mask(badChan) = false;
        vis_artifacts(EEG,oriEEG);
        % vis_artifacts(EEG,oriEEG,'ChannelSubset',1:EEG.nbchan-length(params.heart_channels));
    end

    % Interpolate them
    if EEG.nbchan >30
        EEG = pop_interp(EEG, oriEEG.chanlocs, 'spherical'); % interpolate
        EEG.etc.clean_channel_mask(1:EEG.nbchan) = true;
    else
        warning('Cannot interpolate bad EEG channels reliably with less than 30 channels')
    end

    % % Add ECG channels back?
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
        
        % if params.vis
        %     figure('color','w');
        %     visWeights = 'on';
        % else
        %     visWeights = 'off';
        % end
        % for iChan = 1:EEG.nbchan
        %     [dist(iChan,:),out(iChan,:),rf,w1,w2] = limo_pcout(squeeze(EEG.data(iChan,:,:))','figure',visWeights);
        % end
        % bad_trials = find(isoutlier(sum(~out,1)/EEG.nbchan,'mean'));
        % eegplot(EEG.data(:,:,bad_trials))

        % Detect and remove bad epochs
        badTrials = find_badTrials(EEG,'grubbs', params.vis);
        EEG = pop_rejepoch(EEG, badTrials, 0);

        % % run RMS a 2nd time more conservative in case some were missed
        % badTrials = find_badTrials(EEG,'mean', params.vis);
        % EEG = pop_rejepoch(EEG, badTrials, 0);

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
    
    % Run ICA at effective data rank to control for ghost ICs (Kim et al. 2023). 
    % Use Picard algorithm when possible to increase speed. 
    % use lrate=1e-5 and maxsteps=2000 to obtain reproducible ICA results
    dataRank = sum(eig(cov(double(EEG.data(:,:)'))) > 1E-7);
    if exist('picard.m','file')
        EEG = pop_runica(EEG,'icatype','picard','maxiter',500,'mode','standard', ...
            'pca',dataRank);
    else
        EEG = pop_runica(EEG,'icatype','runica','extended',1, ...
            'pca',dataRank,'lrate',1e-5,'maxsteps',2000);
    end

    % Classify and remove bad components with IClabel
    EEG = pop_iclabel(EEG,'default');
    if strcmp(params.analysis, 'hep') 
        % Do not remove heart components if any
        EEG = pop_icflag(EEG,[NaN NaN; .95 1; .9 1; NaN NaN; .99 1; .99 1; NaN NaN]);
    else
        EEG = pop_icflag(EEG,[NaN NaN; .95 1; .9 1; .99 1; .99 1; .99 1; NaN NaN]);
    end
    badComp = find(EEG.reject.gcompreject);
    EEG = eeg_checkset(EEG);
    
    % Visualize indepent components tagged as bad
    if params.vis
        nComps = size(EEG.icaact,1);
        if nComps >= 20
            pop_selectcomps(EEG,1:20); 
        else
            pop_selectcomps(EEG,1:nComps)
        end
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
