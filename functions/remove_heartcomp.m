%% Remove heart components from EEG signals using ICLabel
%
% Cedric Cannard, 2023

function EEG = remove_heartcomp(EEG, params)

% default ICA algorithm
if ~isfield(params,'icamethod')
    params.icamethod = 1;
end

% default confidence level
if ~isfield(params,'conf_thresh')
    params.conf_thresh = 0.8;  % 80% confidence
else
    if params.conf_thresh > 1  % from GUI is in %
        params.conf_thresh = params.conf_thresh / 100;
    end
end

% Rescale Cardio signal
idx = contains({EEG.chanlocs.labels}, params.heart_channels);
EEG.data(idx,:) = rescale(EEG.data(idx,:), -500, 500);

% smear cardio signal across EEG signals to increase accuracy
if isfield(params,'boost') && params.boost
    % pop_eegplot(EEG,1,1,1);
    disp('Smearing cardiovascular signal across EEG channels to increase classification performance (beta)... ')
    EEG.data = EEG.data - repmat(mean(EEG.data),size(EEG.data,1),1);
    % pop_eegplot(EEG,1,1,1);
end

% Run ICA at effective data rank to control for ghost ICs (Kim et al. 2023). 
% Use Picard algorithm when possible to increase speed. 
% use lrate=1e-5 and maxsteps=2000 to obtain reproducible ICA results
dataRank = sum(eig(cov(double(EEG.data(:,:)'))) > 1E-7);
if params.icamethod == 1
    EEG = pop_runica(EEG,'icatype','picard','maxiter',500,'mode','standard', ...
        'pca',dataRank);
elseif params.icamethod == 2
    EEG = pop_runica(EEG,'icatype','runica','extended',1,'pca',dataRank);
elseif params.icamethod == 3 
    EEG = pop_runica(EEG,'icatype','runica','extended',1, ...
        'pca',dataRank,'lrate',1e-5,'maxsteps',2000);
end

% Classify components with ICLabel
EEG = pop_iclabel(EEG,'default');
EEG = pop_icflag(EEG,[NaN NaN; NaN NaN; NaN NaN; params.conf_thresh 1; NaN NaN; NaN NaN; NaN NaN]); % flag heart components
% pop_selectcomps(EEG,1:EEG.nbchan); colormap('parula');
heart_comp = find(EEG.reject.gcompreject);

% Visualize
if ~isempty(heart_comp)

    fprintf('Number of heart components detected: %g. \n', length(heart_comp));

    % Visualize heart component
    if params.vis_outputs
        if length(heart_comp)==1
            pop_selectcomps(EEG,heart_comp); 
        else
            pop_selectcomps(EEG,1:max(heart_comp)); 
        end
        set(gcf,'Name','Heart component(s) removed','NumberTitle','Off')  % name
        colormap("parula"); pause(0.1)
    end
    
    % Substract heart component from signal
    oriEEG = EEG;
    EEG = pop_subcomp(EEG, heart_comp, 0);
    
    % visualize 
    if params.vis_outputs
        if isfield(EEG.etc, 'clean_channel_mask')
            EEG.etc = rmfield(EEG.etc, 'clean_channel_mask');
        end
        if isfield(oriEEG.etc, 'clean_channel_mask')
            oriEEG.etc = rmfield(oriEEG.etc, 'clean_channel_mask');
        end
        if isfield(EEG.etc, 'clean_sample_mask') 
            EEG.etc = rmfield(EEG.etc, 'clean_sample_mask');
        end
        if isfield(oriEEG.etc, 'clean_sample_mask')
            oriEEG.etc = rmfield(oriEEG.etc, 'clean_sample_mask');
        end
        vis_artifacts(EEG,oriEEG,'ShowSetname',false); 
        set(gcf, 'Toolbar', 'none', 'Menu', 'none','Name', 'Heart components removed', 'NumberTitle', 'Off');                    % remove toolbobar and menu
    end
    
    % Remove CARDIO channel
    % if isfield(params,'keep_heart') && ~params.keep_heart
    EEG = pop_select(EEG,'nochannel', params.heart_channels);
    % end

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
