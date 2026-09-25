% COMPUTE_ASYMMETRY - Alpha asymmetry on all symmetric pairs of electrodes:
% ln(alpha_left + eps) - ln(alpha_right + eps).
%
% eps avoids the log of zero power. Each left electrode is paired with the
% right electrode closest to its mirror position (within 10% of the head
% radius); midline electrodes are excluded. The electrode pairs are
% printed in the command window for visual check.
%
% Usage:
%   [asy, pairLabels, pairNums] = compute_asymmetry(alpha_pwr, norm, chanlocs, vis, tot_pwr)
%
% Inputs:
%   alpha_pwr  - mean alpha power spectral density of each channel (channels x 1).
%                Must be in uV^2/Hz (not log or dB), since the log is taken here.
%   norm       - if true, divide each channel's alpha power by its own total
%                power (tot_pwr) before the log, i.e. relative alpha power
%                (see Allen et al. 2004; Smith et al. 2017)
%   chanlocs   - EEG channel locations (EEGLAB format, with X, Y, Z)
%   vis        - plot asymmetry on a 3D head (left electrodes) (true) or not (false)
%   tot_pwr    - total power of each channel (channels x 1, linear, e.g. mean
%                PSD over the whole frequency range); required if norm is true
%
% Outputs:
%   asy        - asymmetry for each pair (pairs x 1)
%   pairLabels - 'left right' electrode labels of each pair (pairs x 1 cell)
%   pairNums   - [left right] channel indices of each pair (pairs x 2)
%
% Copyright (C) - Cedric Cannard, 2023

function [asy, pairLabels, pairNums] = compute_asymmetry(alpha_pwr, norm, chanlocs, vis, tot_pwr)

if norm && (~exist('tot_pwr','var') || isempty(tot_pwr))
    error("compute_asymmetry: 'tot_pwr' (total power of each channel) is required to normalize asymmetry.")
end

nChan = size(chanlocs,2);
pairNums = nan(nChan,2);
pairLabels = cell(nChan,1);

% Pair each left electrode with the right electrode closest to its mirror
% position (EEGLAB coordinates: X to the nose, Y to the left ear, Z up)
fprintf('Extracting alpha asymmetry on all possible electrode pairs... \n')
try
    XYZ = [[chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]'];
    r = median(sqrt(sum(XYZ.^2,2)), 'omitnan');   % head radius
    tol = 0.02*r;                                 % midline tolerance
    for iChan = find(XYZ(:,2) > tol)'
        d = sqrt(sum((XYZ - XYZ(iChan,:).*[1 -1 1]).^2, 2));
        d(iChan) = Inf;
        [dmin, match] = min(d);
        % mirror electrode must exist (within 10% of the head radius) on the right
        if dmin < 0.1*r && XYZ(match,2) < -tol
            pairNums(iChan,:) = [iChan match];
            pairLabels(iChan,:) = { sprintf('%s %s', chanlocs(iChan).labels, chanlocs(match).labels) };
        end
    end
catch
    warndlg(sprintf("'compute_asymmetry' failed to find the electrode pairs to compute alpha asymmetry. \n\nThis can occur if your dataset contains auxiliary (non-EEG) electrodes, channels without locations, or only one EEG channel."))
end

% Remove empty pairNums
pairNums(cellfun(@isempty,pairLabels),:) = [];
pairLabels(cellfun(@isempty,pairLabels)) = [];

% Remove pairNums with midline electrodes if any made it by mistake
pairNums(contains(pairLabels, 'z'),:) = [];
pairLabels(contains(pairLabels, 'z')) = [];
if size(pairNums,1)~=length(pairLabels)
    warning("Different number of pairs between electroe numbers and labels. There may be an error or NaNs")
end

% Normalize each channel's alpha power by its own total power (relative
% alpha power, see Allen et al. 2004 and Smith et al. 2017)
if norm
    alpha_pwr = alpha_pwr(:) ./ tot_pwr(:);
end

% Natural log of alpha power (+ eps to avoid log of zero)
alpha_pwr = log(alpha_pwr + eps);

nPairs = length(pairLabels);
asy = nan(nPairs,1);
for iPair = 1:nPairs

    % alpha power for each side of the pair
    alpha_left = alpha_pwr(pairNums(iPair,1));
    alpha_right = alpha_pwr(pairNums(iPair,2));

    % Compute asymmetry
    asy(iPair,:) = alpha_left - alpha_right;

end

% 3D plot showing asymmetry on left side of head
if vis
    try
        figure('color','w')
        headplotparams = { 'meshfile','mheadnew.mat','transform',[0.664455 -3.39403 -14.2521 -0.00241453 0.015519 -1.55584 11 10.1455 12],'material','metal' };
        headplot('setup',chanlocs(pairNums(:,1)),'tmp.spl',headplotparams{:}); % Generate temporary spline file
        headplot(asy,'tmp.spl','view',[-85 20],headplotparams{:});  % 3D headplot of asymmetry
        title('Alpha asymmetry')
    catch
        warning("Sorry, 3D headplot failed. Could be because the mesh file was not on the path if using this function outside of BrainBeats.")
    end
end

% Print electrode pairs in command window if some are present
if ~isempty(asy)
    disp('Electrode pairs: ')
    fprintf('   %s \n', pairLabels{:})
    fprintf(['Alpha asymmetry was succesfully computed on %g channel pairs. ' ...
        'Values are ln(left alpha) - ln(right alpha): since alpha is inversely related to cortical activity, positive values reflect greater right-hemispheric activity, and negative values greater left-hemispheric activity \n'], length(asy))
else
    warning("Failed to compute alpha asymmetry on these data")
end


