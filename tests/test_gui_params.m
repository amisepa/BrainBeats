function test_gui_params
% TEST_GUI_PARAMS - Checks of the BrainBeats windows without opening them.
%
% The main window is replaced by preset answers (see BRAINBEATS_MAIN_WINDOW)
% and EEGLAB's inputgui by a stub that returns the default value of every
% control of the parameters window (or the overrides given), and errors if
% the window geometry does not match its controls. For each analysis and
% heart signal, the parameters returned by GETPARAMS_GUI are checked, then
% one HEP analysis is run through the GUI path and its history command is
% run again. Errors on the first failed check.
%
% Usage:
%   test_gui_params          % in MATLAB, EEGLAB on the path
%
% Copyright (C) - Cedric Cannard, 2026

bb = fileparts(fileparts(mfilename('fullpath')));
if ~exist('pop_loadset','file'), eeglab; close; end
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bb));
if ~isempty(old), rmpath(old{:}); end
addpath(bb, fullfile(bb,'functions'));

% inputgui stub, ahead of EEGLAB's on the path for the duration of the test
stub = fullfile(tempdir,'brainbeats_gui_stub');
if ~exist(stub,'dir'), mkdir(stub); end
write_stub(fullfile(stub,'inputgui.m'));
addpath(stub);
cleanup = onCleanup(@() restore(stub));

EEG = pop_loadset('filename','dataset.set','filepath',fullfile(bb,'sample_data'));
base = struct('vis_cleaning',false,'vis_outputs',false,'save',false);

% 1. Defaults of each analysis / heart signal / domain
cases = {
    {'hep','ecg','ECG','channels'}
    {'hep','ecg','ECG','ics'}
    {'hep','ppg','PPG','both'}
    {'features','ecg','ECG',''}
    {'features','ppg','PPG',''}
    {'coherence','ecg','ECG',''}
    };
for c = 1:numel(cases)
    P = run_gui(EEG, base, cases{c});
    assert(strcmp(P.ref,'average') && strcmp(P.clean_method,'asr_ica') && P.linenoise == 50 && strcmp(P.filttype,'noncausal'))
    if strcmp(P.analysis,'hep')
        assert(P.highpass == 0.5 && P.lowpass == 30 && isequal(P.hep_window,[-300 600]) && strcmp(P.hep_baseline,'none'))
        assert(P.hep_surrogates == 100 && strcmp(P.detectMethod,'grubbs') && strcmp(P.heart_removal,'ica') && abs(P.conf_thresh - .75) < 1e-9)
        assert(strcmp(P.hep_level, cases{c}{4}) && isfield(P,'hep_roi') ~= strcmp(P.hep_level,'ics'))
        if strcmp(P.heart_signal,'ppg'), assert(strcmp(P.ppg_transit,'auto') && strcmp(P.ppg_detect_mode,'valleys')); end
    else
        assert(P.highpass == 1 && P.lowpass == 40 && P.asr_cutoff == 30)
    end
    fprintf('OK: %s / %s %s\n', cases{c}{[1 2 4]});
end

% 2. Non-default choices
setappdata(0,'inputgui_override',struct('clean_method',2,'ref',3,'hep_baseline',2,'ppg_transit_mode',2, ...
    'ppg_transit','300','hep_window_mode',2,'linenoise',2,'filttype',2,'detectMethod',1,'ppg_detect_mode',2, ...
    'heart_removal',2));
P = run_gui(EEG, base, {'hep','ppg','PPG','channels'});
assert(strcmp(P.clean_method,'gedai') && strcmp(P.ref,'csd') && strcmp(P.hep_baseline,'regression') && isequal(P.hep_baseline_win,[-150 -50]))
assert(P.ppg_transit == 300 && strcmp(P.hep_window,'adaptive') && P.linenoise == 60 && strcmp(P.filttype,'causal'))
assert(strcmp(P.detectMethod,'quartiles') && strcmp(P.ppg_detect_mode,'peaks') && strcmp(P.heart_removal,'none'))   % PPG: ICA or none
setappdata(0,'inputgui_override',struct('clean_eeg',0));
P = run_gui(EEG, base, {'hep','ecg','ECG','channels'});
assert(~P.clean_eeg && ~isfield(P,'ref') && ~isfield(P,'highpass') && ~isfield(P,'heart_removal'))
setappdata(0,'inputgui_override',struct('heart_removal',2));
P = run_gui(EEG, base, {'hep','ecg','ECG','channels'});
assert(strcmp(P.heart_removal,'ecg_regression'))
fprintf('OK: non-default choices, cleaning off\n');

% 3. A run through the GUI path, and its history command again
setappdata(0,'inputgui_override',struct('clean_eeg',0,'hep_surrogates','10'));
setappdata(0,'brainbeats_main_answers',answers(base, {'hep','ecg','ECG','channels'}));
E = pop_select(EEG,'nochannel',{'PPG'});
[H, com] = brainbeats_process(E);
rmappdata(0,'inputgui_override'); rmappdata(0,'brainbeats_main_answers');
assert(H.trials > 100 && isfield(H.brainbeats,'roi') && isfield(H.brainbeats.roi.tf,'hrpc') && isfield(H.brainbeats,'hrsp'))
EEG = E; %#ok<NASGU> (the command uses EEG)
com = strrep(com, ');', ',''gong'',0);');
eval(com);
assert(EEG.trials == H.trials, 'The history command does not reproduce the GUI run')
fprintf('OK: GUI run (%d epochs) and its history command\n%s\n', H.trials, com);
fprintf('All GUI checks passed.\n')
end

function P = run_gui(EEG, base, c)
setappdata(0,'brainbeats_main_answers',answers(base, c));
[P, abort] = getparams_gui(EEG);
assert(~abort, 'The GUI aborted')
end

function S = answers(base, c)
S = base; S.analysis = c{1}; S.heart_signal = c{2}; S.heart_channels = c(3);
if ~isempty(c{4}), S.hep_level = c{4}; end
end

function restore(stub)
rmpath(stub);
for k = {'inputgui_override' 'brainbeats_main_answers'}
    if isappdata(0,k{1}), rmappdata(0,k{1}); end
end
end

function write_stub(file)
% inputgui stub: default (or overridden) value of each control, by tag
code = {
'function [result, userdat, strhalt, resstruct] = inputgui(geometry, uilist, varargin)'
'nG = sum(cellfun(@numel, geometry));'
'assert(nG == numel(uilist), ''inputgui stub: geometry has %d slots for %d controls'', nG, numel(uilist));'
'result = {}; resstruct = struct(); userdat = []; strhalt = '''';'
'ov = getappdata(0,''inputgui_override'');'
'for i = 1:numel(uilist)'
'    u = uilist{i};'
'    if isempty(u), continue; end'
'    st = u{find(strcmpi(u(1:2:end),''style''))*2};'
'    if ~ismember(st, {''edit'',''popupmenu'',''checkbox'',''listbox''}), continue; end'
'    it = find(strcmpi(u(1:2:end),''tag'')); tag = '''';'
'    if ~isempty(it), tag = u{it*2}; end'
'    if strcmp(st,''edit'')'
'        is = find(strcmpi(u(1:2:end),''string'')); v = '''';'
'        if ~isempty(is), v = u{is*2}; end'
'    else'
'        iv = find(strcmpi(u(1:2:end),''value'')); v = double(strcmp(st,''popupmenu''));'
'        if ~isempty(iv), v = u{iv*2}; end'
'    end'
'    if ~isempty(ov) && isfield(ov, tag), v = ov.(tag); end'
'    result{end+1} = v;'
'    if ~isempty(tag), resstruct.(tag) = v; end'
'end'
};
fid = fopen(file,'w'); fprintf(fid,'%s\n',code{:}); fclose(fid);
rehash
end
