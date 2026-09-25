function EEG = apply_car(EEG, refLabels, chanlocs_template)
% APPLY_CAR - Full-rank common average reference (CAR).
%
% A zero-filled dummy channel (the implicit online reference) is added before
% averaging and removed afterward, so the CAR does not reduce the data rank.
% Optionally, the online reference channel(s) are first added back with
% locations from standard_1005.elc (EEGLAB dipfit) and interpolated.
%
% Usage:
%   EEG = apply_car(EEG)
%   EEG = apply_car(EEG, refLabels, chanlocs_template)
%
% Inputs:
%   EEG               - EEGLAB EEG structure (continuous)
%   refLabels         - (optional) char or cell array of reference label(s) to
%                       recover; labels already in EEG.chanlocs are skipped
%   chanlocs_template - (optional) must be non-empty to enable recovery; only
%                       tested for emptiness (locations come from standard_1005.elc)
% Outputs:
%   EEG - re-referenced EEG structure (recovered reference channels included)
%
% Prints the data rank before and after CAR.
%
% Example:
%   EEG = apply_car(EEG, {'FCz','AFz'}, oriEEG.chanlocs);
%
% References:
%   Makoto's preprocessing pipeline: https://eeglab.ucsd.edu/wiki/Makoto's_preprocessing_pipeline
%   Kim et al. (2023). ICA's bug. Front. Signal Process., 3, 1064138.
%
% Copyright (C) - Cedric Cannard, 2025

rankBefore = sum(eig(cov(double(EEG.data'))) > 1e-7);
fprintf('apply_car: rank before CAR = %d / %d\n', rankBefore, EEG.nbchan);

retainRefs = nargin == 3 && ~isempty(refLabels) && ~isempty(chanlocs_template);

if retainRefs
    if ischar(refLabels) || isstring(refLabels)
        refLabels = cellstr(refLabels);
    end
    refLabels  = refLabels(~ismember(refLabels, {EEG.chanlocs.labels}));
    newChanIdx = zeros(1, numel(refLabels));

    % Add zero-filled channels with placeholder chanlocs
    for i = 1:numel(refLabels)
        EEG.data(end+1, :)       = 0;
        EEG.nbchan               = EEG.nbchan + 1;
        EEG.chanlocs(end+1)      = EEG.chanlocs(1); % temp placeholder, overwritten below
        EEG.chanlocs(end).labels = refLabels{i};
        newChanIdx(i)            = EEG.nbchan;
    end

    % Look up locations in the BEM template; pop_chanedit handles all
    % coordinate conversions (XYZ -> theta/radius/etc.)
    bemPath = which('standard_1005.elc');
    if isempty(bemPath)
        error('apply_car: standard_1005.elc not found. Ensure EEGLAB dipfit plugin is installed.');
    end
    EEG = pop_chanedit(EEG, 'lookup', bemPath);

    % Interpolate the newly added zero-filled channels
    EEG.data = double(EEG.data);
    EEG      = eeg_checkset(EEG);
    EEG      = pop_interp(EEG, newChanIdx, 'spherical');
end

% Add a zero row (the online reference), average-reference, then drop it:
% the remaining channels keep full rank
EEG.data(end+1, :) = 0;
EEG.nbchan         = EEG.nbchan + 1;
EEG.data           = EEG.data - mean(EEG.data, 1);
EEG.data(end, :)   = [];
EEG.nbchan         = EEG.nbchan - 1;

rankAfter = sum(eig(cov(double(EEG.data'))) > 1e-7);
fprintf('apply_car: rank after CAR  = %d / %d\n', rankAfter, EEG.nbchan);

EEG = eeg_checkset(EEG);
end