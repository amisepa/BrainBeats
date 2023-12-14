%% Remove heart components from EEG signals using ICLabel
%
% Cedric Cannard, 2023

function EEG = remove_heartcomp(EEG, CARDIO, params)

% if params.clean_eeg
%     % Filter, re-reference, remove bad channels
%     params.clean_eeg_step = 0;
%     [EEG, params] = clean_eeg(EEG, params); 
% end

% Add CARDIO channels back
EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
for iChan = 1:CARDIO.nbchan
    EEG.chanlocs(end+1).labels = params.heart_channels{iChan};
end
EEG = eeg_checkset(EEG);

% run ICA
dataRank = sum(eig(cov(double(EEG.data(:,:)'))) > 1E-7);
if exist('picard.m','file')
    EEG = pop_runica(EEG,'icatype','picard','maxiter',500,'mode','standard','pca',dataRank);
else
    EEG = pop_runica(EEG,'icatype','runica','extended',1,'pca',dataRank);
end

% Classify components with ICLabel
EEG = pop_iclabel(EEG,'default');
EEG = pop_icflag(EEG,[NaN NaN; NaN NaN; NaN NaN; 0.85 1; NaN NaN; NaN NaN; NaN NaN]); % flag heart components with 85% confidence
% pop_selectcomps(EEG,1:EEG.nbchan); colormap('parula');
heart_comp = find(EEG.reject.gcompreject);

% Visualize
if ~isempty(heart_comp)

    fprintf('Number of heart components detected: %g. \n', length(heart_comp));

    % Visualize heart component
    if params.vis_outputs
        if length(heart_comp)==1
            pop_selectcomps(EEG,heart_comp); colormap("parula")
        else
            pop_selectcomps(EEG,1:max(heart_comp)); colormap("parula")
        end
    end
    
    % Substract heart component from signal
    oriEEG = EEG;
    EEG = pop_subcomp(EEG, heart_comp, 0);
    
    % visualize 
    if params.vis_outputs
        if isfield(EEG.etc, 'clean_channel_mask') %&& length(EEG.etc.clean_channel_mask)<length(oriEEG.etc.clean_channel_mask)
            EEG.etc = rmfield(EEG.etc, 'clean_channel_mask');
            % EEG.etc.clean_channel_mask(end+1:end+CARDIO.nbchan) = true;
        end
        if isfield(oriEEG.etc, 'clean_channel_mask')
            oriEEG.etc = rmfield(oriEEG.etc, 'clean_channel_mask');
        end
        vis_artifacts(EEG,oriEEG,'ShowSetname',false); 
        set(gcf, 'Toolbar', 'none', 'Menu', 'none');                    % remove toolbobar and menu
        set(gcf, 'Name', 'Heart components removed', 'NumberTitle', 'Off')  % name
    end
    
    % Remove CARDIO channel
    if isfield(params,'keep_heart') && ~params.keep_heart
        EEG = pop_select(EEG,'nochannel', params.heart_channels);
    end

else
    fprintf('Sorry, no heart component was detected. Make sure the CARDIO channel you selected is correct. \nYou may try to clean large artifacts in your file to improve ICA performance (or lower the condidence threshold but not recommended). \n')
end

% Store parameters in EEG structure for reporting in publications
EEG.brainbeats.params = params;

% Save
if params.save
    newname = sprintf('%s_no-heart.set', EEG.filename(1:end-4));
    pop_saveset(EEG,'filename',newname,'filepath',EEG.filepath);
end
