% OUTPUT_FILE - File name and folder of a BrainBeats output: next to the input
% file, or, for a dataset that has no file, named after the dataset in the
% current folder.
%
% Usage:
%   [name, folder] = output_file(EEG, suffix)
%
% Inputs:
%   EEG    - EEGLAB dataset
%   suffix - end of the file name, with its extension (e.g. '_HEP.set')
%
% Copyright (C) - Cedric Cannard, 2026

function [name, folder] = output_file(EEG, suffix)
[~, base] = fileparts(char(EEG.filename));
if isempty(base), base = regexprep(char(EEG.setname), '\W+', '_'); end
if isempty(base), base = 'brainbeats'; end
name = [base suffix];
folder = char(EEG.filepath);
if isempty(folder), folder = pwd; end
