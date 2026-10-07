% POP_BRAINBEATS - Open the BrainBeats window (EEGLAB menu 'BrainBeats').
%
% Opens the main window (dataset, analysis, heart signal and channels; see
% BRAINBEATS_MAIN_WINDOW), then the parameters of the analysis, and runs
% BRAINBEATS_PROCESS. A dataset can be loaded from the window if none is
% loaded in EEGLAB.
%
% Usage:
%   [EEG, com] = pop_brainbeats;        % no dataset loaded
%   [EEG, com] = pop_brainbeats(EEG);
%
% Outputs:
%   EEG - processed dataset
%   com - command line reproducing the analysis (EEGLAB history)
%
% Copyright (C) - Cedric Cannard, 2026

function [EEG, com] = pop_brainbeats(EEG)

if nargin < 1, EEG = []; end
[EEG, com] = brainbeats_process(EEG);
