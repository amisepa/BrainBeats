%% Remove heart components from EEG signals using ICLabel
%
% Cedric Cannard, 2023

function EEG = remove_heartcomp(EEG, params)

if params.clean_eeg
    ECG = pop_select(EEG,'channel',params.heart_channels); % export ECG data in separate structure
    EEG = pop_select(EEG,'nochannel',params.heart_channels); % FIXME: remove all non-EEG channels instead

    % Filter, re-reference, remove bad channels
    params.clean_eeg_step = 0;
    [EEG, params] = clean_eeg(EEG, params);

    % Add ECG channels back
    EEG.data(end+1:end+ECG.nbchan,:) = ECG.data;
    EEG.nbchan = EEG.nbchan + ECG.nbchan;
    for iChan = 1:ECG.nbchan
        EEG.chanlocs(end+1).labels = params.heart_channels{iChan};
    end
    EEG = eeg_checkset(EEG);

end

dataRank = sum(eig(cov(double(EEG.data(:,:)'))) > 1E-7);
if exist('picard.m','file')
    EEG = pop_runica(EEG,'icatype','picard','maxiter',500,'mode','standard','pca',dataRank);
else
    EEG = pop_runica(EEG,'icatype','runica','extended',1,'pca',dataRank);
end

EEG = pop_iclabel(EEG,'default');
EEG = pop_icflag(EEG,[NaN NaN; NaN NaN; NaN NaN; 0.95 1; NaN NaN; NaN NaN; NaN NaN]); % flag heart components with 95% confidence
heart_comp = find(EEG.reject.gcompreject);

if ~isempty(heart_comp)

    % Visualize heart component
    if params.vis
        pop_selectcomps(EEG,heart_comp); colormap("parula")
    end

    % Substract heart component from signal
    fprintf('Removing %g heart component(s). \n', length(heart_comp));
    oriEEG = EEG;
    EEG = pop_subcomp(EEG, heart_comp, 0);
    % ADD: option to keep ECG channels by adding them back?

    % visualize with ECG signal to see heartbeats and corresponding
    % contamination in EEG
    if params.vis
        vis_artifacts(EEG,oriEEG);
    end

    % Remove ECG channel
    EEG = pop_select(EEG,'nochannel', params.heart_channels);
else
    fprintf('Sorry, no heart component was detected. Make sure the ECG channel you selected is correct. \nYou may try to clean large artifacts in your file to improve ICA performance (or lower the condidence threshold but not recommended). \n')
end

% Save
if params.save
    newname = sprintf('%s_heart-clean.set', EEG.filename(1:end-4));
    pop_saveset(EEG,'filename',newname,'filepath',EEG.filepath);
end
