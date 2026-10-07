% RUN_GEDAI - Remove artifacts from continuous EEG with GEDAI (generalized
% eigenvalue de-artifacting instrument; Ros et al., 2025).
%
% GEDAI compares the data covariance, in each wavelet band and time window,
% with a leadfield model of brain activity and removes the components that
% do not fit it (eye, muscle, heart and channel artifacts), without ICA. It
% replaces ASR/bad epochs + ICA in CLEAN_EEG when params.clean_method is
% 'gedai'. The output is average-referenced and keeps the input length.
% GEDAI is an external EEGLAB plugin (PolyForm Noncommercial license),
% installed by RUN_CHECKS if missing.
%
% Usage:
%   EEG = run_gedai(EEG, params)
%
% Inputs:
%   EEG    - continuous EEGLAB dataset (EEG channels only, filtered)
%   params - BrainBeats parameters (reads parpool and vis_cleaning)
%
% Output:
%   EEG    - cleaned dataset; GEDAI's scores are in EEG.etc.GEDAI
%
% Notes:
% - Artifact threshold 'auto', 12-cycle windows, 0.5-Hz low cut, and the
%   precomputed 10-05 leadfield (channel labels must be 10-05 names;
%   otherwise the channel positions are used, 'interpolated').
% - Only the first 8 inputs of GEDAI are used: they are the same in GEDAI
%   v1.5 (EEGLAB extension manager) to v1.8.
%
% Reference:
%   Ros, T., et al. (2025). Return of the GEDAI: Unsupervised EEG denoising
%   based on leadfield filtering.
%
% Copyright (C) - Cedric Cannard, 2026

function EEG = run_gedai(EEG, params)

if EEG.trials > 1
    error('run_gedai: GEDAI needs continuous data.')
end
if ~exist('GEDAI','file')
    plugin_askinstall('GEDAI', 'GEDAI', 1);
end
par = isfield(params,'parpool') && ~isempty(params.parpool) && isequal(params.parpool, true);

disp('----------------------------------------------')
fprintf('   Removing artifacts with GEDAI (Ros et al., 2025) \n')
disp('----------------------------------------------')
oriEEG = EEG;
try
    EEG = GEDAI(EEG, 'auto', 12, 0.5, 'precomputed', par, false, Inf);
catch ME
    if contains(ME.identifier, 'ElectrodeLabelsNotFound') || contains(lower(ME.message), 'label')
        warning('GEDAI: channel labels not in its 10-05 template (%s). Using the channel positions instead.', ME.message)
        EEG = GEDAI(oriEEG, 'auto', 12, 0.5, 'interpolated', par, false, Inf);
    else
        rethrow(ME)
    end
end
% GEDAI's own summary: ENOVA = variance of the removed noise as a proportion
% of the original variance (mean across its epochs); SENSAI = denoising
% quality score (%)
if isfield(EEG.etc,'GEDAI')
    g = EEG.etc.GEDAI;
    fprintf('GEDAI: %.0f%% of the variance removed on average (ENOVA), SENSAI score %.0f%%. \n', ...
        100*mean(g.mean_ENOVA), mean(g.SENSAI_score));
end

if isfield(params,'vis_cleaning') && params.vis_cleaning
    vis_artifacts(EEG, oriEEG, 'ShowSetname', false);
    set(gcf,'Toolbar','none','Menu','none','Name','EEG before (red) and after (blue) GEDAI','NumberTitle','Off')
    finish_figure(gcf)
end
