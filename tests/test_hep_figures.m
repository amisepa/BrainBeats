function test_hep_figures
% TEST_HEP_FIGURES - The HEP output figures contain what they should.
%
% Runs a HEP analysis on the sample data (no EEG cleaning, to be fast) with
% the output plots on and the figures invisible, with and without the heart
% channel kept, and checks the content of the figures: the 'HEP' figure must
% hold the image of the single heartbeats (heartbeats x latencies) under the
% HEP of all electrodes, and the 'HEP time-frequency' figure the HRSP and
% HRPC images. erpimage returns without drawing (it only prints a message)
% when it does not accept its inputs, which left the bottom panel empty.
% Errors on the first failed check.
%
% Usage:
%   test_hep_figures          % in MATLAB, EEGLAB on the path
%
% Copyright (C) - Cedric Cannard, 2026

bb = fileparts(fileparts(mfilename('fullpath')));
if ~exist('pop_loadset','file'), eeglab; close; end
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bb));
if ~isempty(old), rmpath(old{:}); end
addpath(bb, fullfile(bb,'functions'));

vis = get(groot,'DefaultFigureVisible');
set(groot,'DefaultFigureVisible','off');
cleanup = onCleanup(@() set(groot,'DefaultFigureVisible',vis));

for keep = [true false]
    close all force
    EEG = pop_loadset('filename','dataset.set','filepath',fullfile(bb,'sample_data'));
    EEG = pop_select(EEG,'nochannel',{'PPG'});
    HEP = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'}, ...
        'clean_eeg',false,'keep_heart',keep,'vis_cleaning',false,'vis_outputs',true, ...
        'save',false,'gong',false);

    f = findall(0,'type','figure','Name','HEP');
    assert(isscalar(f), 'No ''HEP'' figure (keep_heart = %d)', keep)
    im = findall(f,'type','image');
    sz = arrayfun(@(h) size(get(h,'CData'),1:2), im, 'UniformOutput', false);
    ok = cellfun(@(s) s(2) == HEP.pnts && s(1) > 1 && s(1) <= HEP.trials, sz);
    assert(any(ok), 'No image of the single heartbeats in the ''HEP'' figure (keep_heart = %d)', keep)

    f = findall(0,'type','figure','Name','HEP time-frequency');
    assert(isscalar(f) && numel(findall(f,'type','image')) == 2, ...
        'HRSP and HRPC images missing in the ''HEP time-frequency'' figure (keep_heart = %d)', keep)
    fprintf('OK: HEP figures (keep_heart = %d): single heartbeats %d x %d, HRSP, HRPC\n', keep, sz{find(ok,1)});
end
close all force
fprintf('All HEP figure checks passed.\n');
end
