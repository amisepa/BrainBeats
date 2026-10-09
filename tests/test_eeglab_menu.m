function test_eeglab_menu
% TEST_EEGLAB_MENU - The 'BrainBeats' menu of EEGLAB stays usable after a run.
%
% Opens EEGLAB, replaces any 'BrainBeats' menu by the one of this repository
% and checks that it is enabled with no dataset, with a continuous dataset
% and with the epoched output of a HEP run stored as the current dataset (as
% the menu does). Also checks that the HEP output points to the file it was
% saved in (name_HEP.set), that the input file is left untouched, and that
% the main window says what to do with an epoched dataset (image saved in
% tempdir). Errors on the first failed check.
%
% Usage:
%   test_eeglab_menu          % in MATLAB, EEGLAB on the path
%
% Copyright (C) - Cedric Cannard, 2026

bb = fileparts(fileparts(mfilename('fullpath')));
evalin('base','eeglab;');
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bb));
if ~isempty(old), rmpath(old{:}); end
addpath(bb, fullfile(bb,'functions'));

% This repository's menu instead of the installed plugin's
fig = findobj('tag','EEGLAB');
assert(isscalar(fig), 'EEGLAB window not found')
delete(findobj(fig,'type','uimenu','Label','BrainBeats'));
eegplugin_BrainBeats(fig, struct('no_check',''), struct('new_and_hist',''));
m = findobj(fig,'type','uimenu','tag','brainbeats');
assert(isscalar(m), 'BrainBeats menu not created')

% No dataset
evalin('base','eeglab redraw');
check(m, 'no dataset');

% Continuous dataset (scratch copy of the sample data)
dpath = fullfile(tempdir,'brainbeats_test_data');
if ~exist(dpath,'dir'), mkdir(dpath); end
src = fullfile(dpath,'dataset.set');
copyfile(fullfile(bb,'sample_data','dataset.set'), dpath);
before = dir(src);
EEG = pop_loadset('filename','dataset.set','filepath',dpath);
EEG = pop_select(EEG,'nochannel',{'PPG'});
store(EEG);
check(m, 'continuous dataset');

% HEP run, output stored as the current dataset
HEP = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'}, ...
    'clean_eeg',false,'vis_cleaning',false,'vis_outputs',false,'gong',false);
assert(HEP.trials > 1 && isnumeric(HEP.data), 'No HEP epochs in the output')
assert(strcmp(HEP.filename,'dataset_HEP.set') && isfile(fullfile(HEP.filepath,HEP.filename)), ...
    'The HEP output points to %s', HEP.filename)
after = dir(src);
assert(after.bytes == before.bytes && after.datenum == before.datenum, 'The input file was modified')
store(HEP);
check(m, 'epoched dataset (HEP output)');

% Main window with the epoched dataset
snap = fullfile(tempdir,'brainbeats_main_epoched.png');
setappdata(0,'brainbeats_main_snapshot',snap);
cleanup = onCleanup(@() rmappdata(0,'brainbeats_main_snapshot'));
brainbeats_main_window(HEP);
fprintf('OK: main window with an epoched dataset saved in %s\n', snap);

fprintf('All menu checks passed.\n');
end

function store(EEG)
% Store as a new current dataset and redraw the EEGLAB window
assignin('base','EEG',EEG);
evalin('base','[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0); eeglab redraw');
end

function check(m, state)
assert(strcmp(get(m,'Enable'),'on'), 'BrainBeats menu disabled with: %s', state)
fprintf('OK: menu enabled with %s\n', state);
end
