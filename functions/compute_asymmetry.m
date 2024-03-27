%% Compute alpha asymmetry on all symmetric pairs of electrodes, using
% log(alpha_pwr_left + constant) - log(alpha_pwr_right + constant). The
% constant (eps) is to avoid issues when computing log of zero or negative 
% power values. The labels of the pairNums of electrodes are printed in 
% command window for visual check, and a 3D head plot displays the asymmetry 
% output (left side only since it's a relative difference).
%
% INPUTS:
%   alpha_pwr   - mean aplha power spectral density (PSD) for each EEG 
%               channel (channel x power). Warning: Must be in μV^2/Hz to avoid log
%               of already log power. 
%   norm        - normalize (1) or not (0) by dividing channel alpha power 
%               by total alpha power
%   chanlocs    - EEG electrode locations in the EEGLAB format
%   vis         - visualize the outputs on a 3D head image (1) or not (0)
% 
% OUTPUTS:
%   asy         - asymmetry values for each pair
%   pairLabels  - electrode labels of each pair
%   pairNums    - electrode number of each pair
% 
% EXAMPLE:
%   [asy, pairLabels, pairNums] = compute_asymmetry(alpha_pwr, norm, chanlocs, vis)
% 
% Copyright (C) - BrainBeats - Cedric Cannard - 2023

function [asy, pairLabels, pairNums] = compute_asymmetry(alpha_pwr, norm, chanlocs, vis)

nChan = size(chanlocs,2);
elec_dist = nan(nChan,1);
pairNums = nan(nChan,2);
pairLabels = cell(nChan,1);

% Calculate theta from XYZ if not present
ThetaRad = cart2sph([chanlocs.X]',[chanlocs.Y]',[chanlocs.Z]');
theta = -ThetaRad * 180 / pi;

fprintf('Extracting (z-normalized) alpha asymmetry on all possible pairNums... \n')
try
    for iChan = 1:nChan
    
        % if not a left electrode, skip
        if theta(iChan) >= 0, continue; end
    
        % Find matching right electrode using X distance
        for iChan2 = 1:nChan
            if iChan2 ~= iChan
                elec_dist(iChan2,:) = abs( abs(theta(iChan)) - abs(theta(iChan2)) );
            end
        end
        [~, match] = mink(elec_dist,2);
    
        % Fix: remove match that have distance greater than 1, useful for 
        % low-density montages or to prevent errors below with labels like TP9-TP10
        match(elec_dist(match) > 1) = [];
    
        % Fix: sometimes get the wrong one depending on montage where some
        % other electrodes have shorter distance
        match_labels = { chanlocs(match).labels };
        if length(match_labels) > 1
            idx = length(chanlocs(iChan).labels) == cellfun('length',match_labels)';
            if sum(idx)==1
                match = match(idx);
            else
                match = match(1);
            end
        end
    
        % Ensure pairNums are always left-right (i.e. 2nd column should always
        % be right hemisphere) and store
        % if rem(str2double(pairLabels{iChan}(end)),2) ~= 0 % only works with 10-20 labels
        if theta(match) >= 0
            pairNums(iChan,:) = [iChan match];
            pairLabels(iChan,:) = { sprintf('%s %s', chanlocs(iChan).labels, chanlocs(match).labels) };
        end
    end
catch
    errordlg(sprintf("'get_eef_features.m' failed to find all channel pairs to compute alpha asymmetry. \n\nThis is likely because your dataset still contains a non-EEG electrode with a label other than ECG, PPG, or AUX (which would have been detected). \n\nPlease inspect your electrode labels and remove any electrode that should not be there prior to launching BrainBeats."))
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

% Remove errors/duplicate
% pairLabels(pairNums(:,1) == 0) = [];
% pairNums(pairNums(:,1) == 0,:) = [];
% [pairLabels, idx] = unique(pairLabels);
% pairNums = pairNums(idx,:);

% natural logarithm alpha power + eps constant to avoid issues with zero or 
% negative values
alpha_pwr = log(alpha_pwr + eps);           % standard formula
% alpha_pwr = 10*log10(alpha_pwr + eps);      % decibels

nPairs = length(pairLabels);
asy = nan(nPairs,1);
for iPair = 1:nPairs
    
    % alpha power for each side of the pair
    alpha_left = alpha_pwr(pairNums(iPair,1));
    alpha_right = alpha_pwr(pairNums(iPair,2));

    % Normalize by correcting for overall alpha power (see Allen et al. 2004 and Smith et al. 2017)
    if norm
        alpha_left = alpha_left / sum(mean(alpha_pwr,2));
        alpha_right = alpha_right / sum(mean(alpha_pwr,2));
    end

    % Compute asymmetry
    asy(iPair,:) = alpha_left - alpha_right;

end

% 3D plot showing asymmetry on left side of head
if vis
    try
        figure('color','w')
        headplotparams = { 'meshfile','mheadnew.mat','transform',[0.664455 -3.39403 -14.2521 -0.00241453 0.015519 -1.55584 11 10.1455 12],'material','metal' };
        % headplotparams = {'meshfile','colin27headmesh.mat','transform',[0 -13 0 0.1 0 -1.57 11.7 12.5 12],'material','metal' };
        headplot('setup',chanlocs(pairNums(:,1)),'tmp.spl',headplotparams{:}); % Generate temporary spline file
        headplot(asy,'tmp.spl','view',[-85 20],headplotparams{:});  % 3D headplot of asymmetry
        % headplot('setup',chanlocs,'tmp.spl',headplotparams{:}); % Generate temporary spline file
        % headplot(log(alpha_pwr),'tmp.spl','view',[-85 20],headplotparams{:});  % 3D headplot of asymmetry
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
        'Positive values from the output reflect greater left-hemispheric activity, and negative values reflect greater right-hemispheric activity \n'], length(asy))
else
    warning("Failed to compute alpha asymmetry on these data")
end


