% REMOVE_HEARTCOMP - Remove heart components from EEG signals with ICA and ICLabel.
%
% The heart channel(s) are rescaled to +/-500 and kept in the ICA so that
% ICLabel can better isolate the cardiac field artifact (CFA). ICA is run at
% the effective data rank (Kim et al. 2023), ICs labeled heart with
% confidence >= params.conf_thresh are subtracted, and the heart channel(s)
% are then removed unless params.keep_heart (kept in rescaled units).
%
% The CFA left in the EEG is reported before and after: RMS over EEG
% channels of the heartbeat-locked average from -50 to 100 ms (baseline -200
% to -100 ms), with the heartbeats detected on the heart channel. It also
% contains genuine heartbeat-evoked activity, so it does not reach zero.
%
% Usage:
%   EEG = remove_heartcomp(EEG, params)
%
% Inputs:
%   EEG    - EEGLAB EEG structure with EEG and heart channels (continuous)
%   params - BrainBeats parameters:
%            heart_channels - label(s) of the heart channel(s)
%            conf_thresh    - ICLabel heart confidence, 0-1 or in % (default .9)
%            icamethod      - 1 = Picard, 2 = extended Infomax (run_checks
%                             default), 3 = replicable Infomax (lrate 1e-5,
%                             maxsteps 2000); 1 if the field is missing
%            keep_heart     - true to keep the heart channel(s) (default false)
%            vis_outputs, save
% Outputs:
%   EEG    - cleaned EEG, with EEG.brainbeats.parameters and, in
%            EEG.brainbeats.preprocessings, removed_heart_components (IC
%            indices) and cfa_before/cfa_after (uV).
%            Saved as <filename>_no-heart.set if params.save.
%
% Copyright (C) - Cedric Cannard, 2023

function EEG = remove_heartcomp(EEG, params)

% default ICA algorithm
if ~isfield(params,'icamethod')
    params.icamethod = 1;
end

% default confidence level
if ~isfield(params,'conf_thresh')
    params.conf_thresh = 0.9;  % 90% confidence
else
    if params.conf_thresh > 1  % from GUI is in %
        params.conf_thresh = params.conf_thresh / 100;
    end
end

% Rescale cardio signal to an EEG-like amplitude range (units vary across devices)
idx = ismember(lower({EEG.chanlocs.labels}), lower(params.heart_channels));  % exact labels, not substrings
EEG.data(idx,:) = rescale(EEG.data(idx,:), -500, 500);

% Heartbeats on the (first) heart channel, for the CFA readout
hIdx = find(idx, 1);
beatParams = params;
beatParams.fs = EEG.srate;
try
    [~, ~, beats] = get_RR(EEG.data(hIdx,:), EEG.times, beatParams);
    cfaBefore = cfa_amplitude(EEG.data(~idx,:), beats, EEG.srate);
catch
    beats = []; cfaBefore = NaN;
end

% Run ICA at the effective data rank to avoid ghost ICs (Kim et al. 2023).
% Method 3 uses lrate = 1e-5 and maxsteps = 2000 for reproducible results.
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

% Classify components with ICLabel and flag heart components only
% (rows: brain, muscle, eye, heart, line noise, channel noise, other)
EEG = pop_iclabel(EEG,'default');
EEG = pop_icflag(EEG,[NaN NaN; NaN NaN; NaN NaN; params.conf_thresh 1; NaN NaN; NaN NaN; NaN NaN]); % flag heart components
heart_comp = find(EEG.reject.gcompreject);

if ~isempty(heart_comp)

    fprintf('Number of heart components detected: %g. \n', length(heart_comp));

    % Plot heart component(s) (all ICs up to the last heart one if several)
    if params.vis_outputs
        if length(heart_comp)==1
            pop_selectcomps(EEG,heart_comp); 
        else
            pop_selectcomps(EEG,1:max(heart_comp)); 
        end
        set(gcf,'Name','Heart component(s) removed','NumberTitle','Off')  % name
        colormap('parula'); pause(0.1)
    end
    
    % Subtract heart component(s) from the signals
    oriEEG = EEG;
    EEG = pop_subcomp(EEG, heart_comp, 0);

    % Plot data before/after (clean_rawdata masks removed first so that
    % vis_artifacts compares the two datasets directly)
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
        set(gcf, 'Toolbar', 'none', 'Menu', 'none','Name', 'Heart components removed', 'NumberTitle', 'Off');                    % remove toolbar and menu
    end
    
    cfaAfter = cfa_amplitude(EEG.data(~idx,:), beats, EEG.srate);
    fprintf('Heart-locked EEG amplitude (-50 to 100 ms): %.2f uV before, %.2f uV after (%.0f%% reduction). \n', ...
        cfaBefore, cfaAfter, 100*(1 - cfaAfter/cfaBefore));
else
    cfaAfter = cfaBefore;
    fprintf('Sorry, no heart component was detected. Make sure the CARDIO channel you selected is correct. \nYou may try to clean large artifacts in your file to improve ICA performance (or lower the condidence threshold but not recommended). \n')
end

% Remove CARDIO channel (whether or not heart components were found)
if ~isfield(params,'keep_heart') || ~params.keep_heart
    EEG = pop_select(EEG,'nochannel', params.heart_channels);
end

% Store parameters in EEG structure for reporting in publications
EEG.brainbeats.parameters = params;
EEG.brainbeats.preprocessings.removed_heart_components = heart_comp;
EEG.brainbeats.preprocessings.cfa_before = cfaBefore;
EEG.brainbeats.preprocessings.cfa_after = cfaAfter;

% Save
if params.save
    newname = sprintf('%s_no-heart.set', EEG.filename(1:end-4));
    pop_saveset(EEG,'filename',newname,'filepath',EEG.filepath);
end


% RMS over channels of the heartbeat-locked average from -50 to 100 ms
% (baseline -200 to -100 ms), in the units of the data
function a = cfa_amplitude(D, beats, fs)
w = round(-0.2*fs):round(0.4*fs);
beats = beats(beats + w(1) > 0 & beats + w(end) <= size(D,2));
if numel(beats) < 10
    a = NaN; return
end
avg = zeros(size(D,1), numel(w));
for i = 1:numel(beats)
    avg = avg + double(D(:, beats(i) + w));
end
avg = avg / numel(beats);
t = w / fs;
avg = avg - mean(avg(:, t < -0.1), 2);
a = sqrt(mean(mean(avg(:, t >= -0.05 & t <= 0.1).^2)));
