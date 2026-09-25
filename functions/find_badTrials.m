% FIND_BADTRIALS - Detect bad epochs from their amplitude and high-frequency content.
%
% Usage:
%   badTrials = find_badTrials(EEG, method, vis)
%
% For each epoch, two metrics are computed: amplitude (RMS across channels of
% each channel's RMS) and high-frequency noise (RMS across channels of each
% channel's mean absolute deviation over time of the data minus its 45-Hz
% low-passed version). Epochs that are outliers on either metric
% (isoutlier) are flagged.
%
% Inputs:
%   EEG       - epoched EEGLAB EEG structure
%   method    - isoutlier method: 'median' (more aggressive), 'quartiles',
%               'grubbs' (moderate, clean_eeg default), 'mean' (more lax)
%   vis       - true to plot the flagged epochs (eegplot)
% Outputs:
%   badTrials - indices of bad epochs (row vector)
%
% Copyright (C) - Cedric Cannard, 2022

function badTrials = find_badTrials(EEG,method,vis)

disp('Detecting bad trials...')
% Low-pass FIR (order 100, pass 0-45 Hz, stop from 50 Hz)
b = design_fir(100,[2*[0 45 50]/EEG.srate 1],[1 1 0 0]);
sig_amp = nan(1,size(EEG.data,3));
sig_snr = nan(1,size(EEG.data,3));
for iEpoch = 1:size(EEG.data,3)
    sig_amp(:,iEpoch) = rms(rms(squeeze(EEG.data(:,:,iEpoch)),2));
    tmp = filtfilt_fast(b,1, squeeze(EEG.data(:,:,iEpoch))');
    sig_snr(:,iEpoch) = rms(mad(squeeze(EEG.data(:,:,iEpoch)) - tmp', 0, 2));   % deviation over time (dim 2) for each channel
end
badRMS = isoutlier(sig_amp,method);
badSNR = isoutlier(sig_snr,method);
badTrials = unique([find(badRMS) find(badSNR)]);

if ~isempty(badTrials) 
    if vis
        % bad epochs only, with their own events; eegplot sets the channel
        % spacing from the data (bad epochs are larger than usual)
        BAD = pop_select(EEG, 'trial', badTrials);
        eegplot(BAD.data,'srate',BAD.srate,'events',BAD.event,'eloc_file',BAD.chanlocs, ...
            'winlength',min(5,numel(badTrials)),'title','Epochs removed','plottitle','Bad epochs');
        set(gcf,'Menu', 'none','Name','Bad epochs removed','NumberTitle','Off')
        finish_figure(gcf)
    end
else
    disp("No bad trials detected")
end

fprintf('Trials detected: %g \n', length(badTrials));
