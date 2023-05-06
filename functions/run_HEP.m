% Run heartbeat evoked potentials (HEP)
%
% Cedric Cannard, 2023

function run_HEP(EEG, params, Rpeaks)

% Add heartbeat markers in the signals, taking into account
% existing events
nEv = length(EEG.event);
urevents = num2cell(nEv+1:nEv+length(Rpeaks));
evt = num2cell(Rpeaks);
types = repmat({'R-peak'},1,length(evt));
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).latency] = evt{:};       % assign latencies
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).type] = types{:};        % assign types
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).urevent] = urevents{:};  % assign event index
EEG = eeg_checkset(EEG);

% Calculate minimum distance between R peaks (in s)
minEpoch = min(diff(Rpeaks)/EEG.srate);

% Epoch
EEG = pop_epoch(EEG,{},[-minEpoch minEpoch],'epochinfo','yes');

if params.clean_eeg

    % Remove bad trials
    disp('Looking for bad trials...')
    b = design_fir(100,[2*[0 45 50]/EEG.srate 1],[1 1 0 0]);
    sigRMS = rms(squeeze(rms(EEG.data,2)),1);
    badRMS = isoutlier(sigRMS,'mean');
    snr = nan(1,size(EEG.data,3));
    for iEpoch = 1:size(EEG.data,3)
        tmp = filtfilt_fast(b,1, squeeze(EEG.data(:,:,iEpoch))');
        snr(:,iEpoch) = rms(mad(squeeze(EEG.data(:,:,iEpoch)) - tmp'));
    end
    badSNR = isoutlier(snr,'grubbs');
    badTrials = unique([find(badRMS) find(badSNR)]);
    % pop_eegplot(EEG,1,1,1);
    EEG = pop_rejepoch(EEG, badTrials, 0);

    % Run ICAlabel
    dataRank = sum(eig(cov(double(EEG.data(:,:)'))) > 1E-7);
    if exist('picard.m','file')
        EEG = pop_runica(EEG,'icatype','picard','maxiter',500,'mode','standard','pca',dataRank);
    else
        EEG = pop_runica(EEG,'icatype','runica','extended',1,'pca',dataRank);
    end

    EEG = pop_iclabel(EEG,'default');
    EEG = pop_icflag(EEG,[NaN NaN; .95 1; .95 1; NaN NaN; NaN NaN; NaN NaN; NaN NaN]);
    badComp = find(EEG.reject.gcompreject);
    EEG = eeg_checkset(EEG);
    if params.vis, pop_selectcomps(EEG,1:6); end
    if ~isempty(badComp)
        fprintf('Removing %g bad component(s). \n', length(badComp));
        oriEEG = EEG;
        EEG = pop_subcomp(EEG, badComp, 0);
        % if params.vis, vis_artifacts(EEG,oriEEG); end
    end
end

if params.vis
    pop_eegplot(EEG,1,1,1);
    % eegplot(EEG.data,'winlength',15,'srate',EEG.srate,'events',EEG.event,'spacing',100);
end

% Save
newname = sprintf('%s_HEP.set', EEG.filename(1:end-4));
pop_saveset(EEG,'filename',newname,'filepath',EEG.filepath);

