function EEG = apply_car(EEG, refLabel)
% apply_car - Common average reference with rank preservation
%
%   This method adds a zero-filled surrogate reference channel,
%   computes the CAR including that channel, then removes it.
%
%   EEG = apply_car(EEG, refLabel)
%
% Inputs:
%   EEG      - EEGLAB EEG struct
%   refLabel - Label for the surrogate reference channel (string)
%
% Example:
%   EEG = apply_car(EEG, 'initialReference');
% 
% Reference: 
%   Kim, H., Luo, J., Chu, S., Cannard, C., Hoffmann, S., & Miyakoshi, M. (2023). 
%   ICA's bug: How ghost ICs emerge from effective rank deficiency caused by EEG 
%   electrode interpolation and incorrect re-referencing. Frontiers in Signal Processing, 3, 1064138.


disp("Performing (modified) common average reference (CAR; preserving full data rank by adding a zero-filled channel). For more detail, see: https://www.frontiersin.org/journals/signal-processing/articles/10.3389/frsip.2023.1064138/full");

if nargin < 2
    refLabel = 'initialReference';
end

% 1. Append zero-filled surrogate channel
EEG.data(end+1, :) = 0; % zero signal
EEG.nbchan = EEG.nbchan + 1;

% Clone an existing chanlocs entry for metadata
dummyLoc = EEG.chanlocs(1);
dummyLoc.labels = refLabel;
EEG.chanlocs(end+1) = dummyLoc;

% 2. Apply CAR manually
avgSignal = mean(EEG.data, 1);
EEG.data = EEG.data - avgSignal;

% 3. Remove surrogate channel
keepIdx = find(~strcmp({EEG.chanlocs.labels}, refLabel));
EEG.data = EEG.data(keepIdx, :);
EEG.chanlocs = EEG.chanlocs(keepIdx);
EEG.nbchan = numel(EEG.chanlocs);

