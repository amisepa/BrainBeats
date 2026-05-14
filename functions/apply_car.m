function EEG = apply_car(EEG, refLabels, chanlocs_template)
% apply_car - Full-rank common average reference (CAR)
%
%   Optionally recovers online reference channel(s) by adding them with proper
%   locations and interpolating them before CAR. A zero-filled dummy channel
%   is always added for rank preservation during CAR and removed afterward.
%
%   EEG = apply_car(EEG)
%   EEG = apply_car(EEG, refLabels, chanlocs_template)
%
% Inputs:
%   EEG               - EEGLAB EEG struct
%   refLabels         - string or cell array of reference label(s) to recover
%   chanlocs_template - chanlocs struct or path to electrode file (used as
%                       fallback; BEM standard_1005.elc is primary source).
%
% Outputs:
%   EEG - re-referenced EEG struct. Reference channel(s) are interpolated
%         and retained before CAR is applied.
%
% Usage:
%   EEG = apply_car(EEG);
%   EEG = apply_car(EEG, {'FCz','AFz'}, oriEEG.chanlocs);
%
% References:
%   Makoto's preprocessing pipeline: https://eeglab.ucsd.edu/wiki/Makoto's_preprocessing_pipeline
%   Kim et al. (2023). ICA's bug. Front. Signal Process., 3, 1064138.

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

    % Look up correct locations using pop_chanedit with BEM template
    % This handles all coordinate conversions (XYZ -> theta/radius/etc.)
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

% Always add dummy zero row for rank preservation, apply CAR, remove dummy
EEG.data(end+1, :) = 0;
EEG.nbchan         = EEG.nbchan + 1;
EEG.data           = EEG.data - mean(EEG.data, 1);
EEG.data(end, :)   = [];
EEG.nbchan         = EEG.nbchan - 1;

rankAfter = sum(eig(cov(double(EEG.data'))) > 1e-7);
fprintf('apply_car: rank after CAR  = %d / %d\n', rankAfter, EEG.nbchan);

EEG = eeg_checkset(EEG);
end