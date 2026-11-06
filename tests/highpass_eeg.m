function EEG = highpass_eeg(EEG, cutoff)
% HIGHPASS_EEG - High-pass filter the EEG channels of a test dataset, heart
% channels (ECG, PPG) left as they are. BrainBeats refuses unfiltered EEG for
% HEP without 'clean_eeg', since the HEP is not baseline-corrected.
%
% Usage:
%   EEG = highpass_eeg(EEG)          % 0.5 Hz
%   EEG = highpass_eeg(EEG, cutoff)
%
% Copyright (C) - Cedric Cannard, 2026

if nargin < 2, cutoff = 0.5; end
F = pop_eegfiltnew(EEG, 'locutoff', cutoff);
iE = ~contains(lower({EEG.chanlocs.labels}), {'ecg' 'ekg' 'ppg'});
EEG.data(iE,:) = F.data(iE,:);
end
