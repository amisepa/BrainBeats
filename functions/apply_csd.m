% APPLY_CSD - Surface Laplacian (current source density) of the EEG channels
% of a dataset, after artifact removal ('ref','csd').
%
% The EEG is average-referenced during cleaning (ICA and ICLabel need
% average-referenced data), and the surface Laplacian is applied here, at
% the end, to the EEG channels only (heart channels are left as they are).
% If the transform fails (e.g. channel labels that are not 10-05 names), a
% warning is given and the data stay average-referenced.
%
% Usage:
%   [EEG, C] = apply_csd(EEG)
%   [EEG, C] = apply_csd(EEG, exclude)
%
% Inputs:
%   EEG     - EEGLAB dataset (continuous or epoched)
%   exclude - labels of channels to leave untouched (e.g. heart channels)
%
% Outputs:
%   EEG - with the EEG channels CSD-transformed (EEG.ref = 'csd-transform')
%   C   - transform matrix (EEG channels x EEG channels), [] if it failed
%
% See also: CSD_TRANSFORM
%
% Copyright (C) - Cedric Cannard, 2026

function [EEG, C] = apply_csd(EEG, exclude)

if nargin < 2, exclude = {}; end
idx = find(~ismember(lower({EEG.chanlocs.labels}), lower(cellstr(exclude))));
C = [];
try
    sub = pop_select(EEG, 'channel', idx);
    [~, C] = csd_transform(sub);
    sz = size(EEG.data);
    X = reshape(double(EEG.data), sz(1), []);
    X(idx,:) = C * X(idx,:);
    EEG.data = reshape(X, sz);
    EEG.ref = 'csd-transform';
    if isfield(EEG.chanlocs,'ref')   % eeg_checkset copies chanlocs(1).ref into EEG.ref
        [EEG.chanlocs(idx).ref] = deal('csd-transform');
    end
    EEG.icaact = [];
    fprintf('Surface Laplacian (current source density) applied to %d EEG channels. \n', numel(idx));
catch ME
    warning('The surface Laplacian (CSD) failed (%s): the data stay average-referenced.', ME.message)
end
