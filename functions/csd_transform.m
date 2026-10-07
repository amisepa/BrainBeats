% CSD_TRANSFORM - Current source density (CSD; surface Laplacian) of EEG data,
% with spherical splines (Perrin et al., 1989; Kayser & Tenke, 2006).
%
% Reference-free: the CSD estimates the radial current at each electrode from
% the curvature of the potential field, which sharpens the topographies and
% reduces volume conduction (e.g. of the cardiac field artifact) between
% neighboring electrodes. Channel positions are taken from the 10-05
% template (standard_1005.ced) by label, so every channel must have a
% standard 10-05 name.
%
% Usage:
%   EEG = csd_transform(EEG)
%   [EEG, C] = csd_transform(EEG, chanlocfile, mcont, smoothl, headrad)
%
% Inputs:
%   EEG         - EEGLAB dataset (continuous or epoched), EEG channels only
%   chanlocfile - electrode file (.ced, .xyz, .locs or .csd; default
%                 standard_1005.ced next to this file)
%   mcont       - spline flexibility m, 2-10 (default 4)
%   smoothl     - smoothing constant lambda (default 1e-5)
%   headrad     - head radius in cm (default 10: uV/cm^2)
%
% Outputs:
%   EEG - CSD-transformed dataset (EEG.ref = 'csd-transform')
%   C   - the transform as a channels x channels matrix (the CSD is linear:
%         csd = C * data), to apply the same transform to other data of the
%         same montage
%
% From csd_transfrom (github.com/amisepa/csd_transfrom, GPL-3.0), which uses
% GetGH and CSD by Juergen Kayser (CSD toolbox, GPL).
%
% Copyright (C) - Cedric Cannard, 2025

function [EEG, C] = csd_transform(EEG, chanlocfile, mcont, smoothl, headrad)

if ~exist('mcont','var') || isempty(mcont), mcont = 4; end
if ~exist('smoothl','var') || isempty(smoothl), smoothl = 1.0e-5; end
if ~exist('headrad','var') || isempty(headrad), headrad = 10; end
if isfield(EEG,'ref') && ischar(EEG.ref) && strcmpi(EEG.ref,'csd-transform')
    warning('These data are already CSD-transformed.');
end

% Electrode positions (theta, phi) from the location file, in EEG order
if ~exist('chanlocfile','var') || isempty(chanlocfile)
    chanlocfile = fullfile(fileparts(which('csd_transform.m')), 'standard_1005.ced');
end
if ~exist(chanlocfile,'file')
    error('Channel location file not found: %s', chanlocfile);
end
[~,~,ext] = fileparts(chanlocfile);
[labels, theta, phi] = read_chanlocs(chanlocfile, lower(strrep(ext,'.','')));
eegLabels = {EEG.chanlocs.labels};
[found, loc] = ismember(lower(eegLabels), lower(labels(:)));
if ~all(found)
    error('csd_transform: channel labels not in %s: %s', chanlocfile, strjoin(eegLabels(~found), ', '));
end
M.lab   = labels(loc);
M.theta = theta(loc);
M.phi   = phi(loc);

% G and H matrices, then the linear transform (identity input = the matrix)
[G, H] = GetGH(M, mcont);
nChan = numel(eegLabels);
C = current_source_density(eye(nChan), G, H, smoothl, headrad);

% Apply to the data (continuous or epoched)
sz = size(EEG.data);
EEG.data = reshape(C * double(reshape(EEG.data, nChan, [])), sz);
if isfield(EEG,'chaninfo') && isfield(EEG.chaninfo,'nosedir')
    EEG.chaninfo = struct('nosedir', EEG.chaninfo.nosedir);
end
EEG.ref = 'csd-transform';
EEG.icaact = [];
EEG = eeg_checkset(EEG);
