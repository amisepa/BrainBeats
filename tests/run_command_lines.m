function results = run_command_lines(ids)
% RUN_COMMAND_LINES - Run the documented command lines as they are written.
%
% Runs each section of brainbeats_tutorial.m (text read from the file, so
% the commands are the ones users copy), then the examples of
% 'help brainbeats_process' and README.md, and prints PASS/FAIL for each.
% run_tutorial_headless has its own copy of the calls, which can drift from
% the tutorial; this test cannot.
%
% Usage:
%   results = run_command_lines;        % in MATLAB, EEGLAB on the path
%   results = run_command_lines(1:3);   % a subset of the sections
%
% Input:
%   ids     - indices of the sections to run (default: all)
%
% Output:
%   results - struct array with fields name, status ('PASS'/'FAIL'),
%             seconds and message (error and stack)
%
% Notes:
% - Plots and the end gong are forced off by a brainbeats_process wrapper
%   placed ahead on the path (drawing can hang under 'matlab -batch'), and
%   errordlg/warndlg are replaced by stubs that print the message.
% - The sample data are read from a copy in tempdir, so 'save' never writes
%   into the repo. The first tutorial section (EEGLAB start, cd) and the
%   last one (BrainBeats window) are skipped.
%
% Copyright (C) - Cedric Cannard, 2026

% Paths: EEGLAB, then this repository instead of any other BrainBeats copy
bbroot = fileparts(fileparts(mfilename('fullpath')));
if ~exist('pop_loadset','file'), eeglab; close; end
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bbroot));
if ~isempty(old), rmpath(old{:}); end
addpath(bbroot, fullfile(bbroot,'functions'));
fprintf('Testing %s\n', which('brainbeats_process'));

% Scratch copy of the sample data
dpath = fullfile(tempdir,'brainbeats_test_data');
if ~exist(dpath,'dir'), mkdir(dpath); end
copyfile(fullfile(bbroot,'sample_data','dataset.set'), dpath);

% Sections of the tutorial that call brainbeats_process, reading the data
% from the scratch copy
txt = fileread(fullfile(bbroot,'brainbeats_tutorial.m'));
sec = regexp(txt, '(?m)^%%', 'split');
sec = sec(contains(sec,'brainbeats_process(EEG'));
S = struct('name',{},'code',{});
for k = 1:numel(sec)
    name = strtrim(regexp(sec{k}, '^[^\r\n]*', 'match', 'once'));
    S(end+1) = struct('name',['tutorial: ' name], 'code', ...
        strrep(['%' sec{k}], 'fullfile(main_path,''sample_data'')', 'dpath')); %#ok<AGROW>
end

% Examples of the brainbeats_process help (the README one is the first)
load_ecg = 'EEG = pop_loadset(''filename'',''dataset.set'',''filepath'',dpath); ';
S(end+1) = struct('name','help/README: HEP, ECG, 100 surrogates','code', ...
    [load_ecg 'EEG = brainbeats_process(EEG, ''analysis'',''hep'', ''heart_signal'',''ecg'', ' ...
    '''heart_channels'',{''ECG''}, ''clean_eeg'',true, ''hep_surrogates'',100);']);
S(end+1) = struct('name','help: HEP, PPG, channels and ICs','code', ...
    [load_ecg 'EEG = brainbeats_process(EEG, ''analysis'',''hep'', ''heart_signal'',''ppg'', ' ...
    '''heart_channels'',{''PPG''}, ''clean_eeg'',true, ''hep_level'',''both'');']);
helptxt = help('brainbeats_process');
assert(contains(helptxt,'''hep_surrogates'',100);') && contains(helptxt,'''hep_level'',''both'');'), ...
    'The examples of the brainbeats_process help changed: update this test')

% Stubs: dialogs printed, brainbeats_process with plots and gong off
stubs = fullfile(tempdir,'brainbeats_cmdline_stubs');
if ~exist(stubs,'dir'), mkdir(stubs); end
setappdata(0,'brainbeats_real', @brainbeats_process);
writestub(fullfile(stubs,'errordlg.m'), 'errordlg');
writestub(fullfile(stubs,'warndlg.m'),  'warndlg');
writewrapper(fullfile(stubs,'brainbeats_process.m'));
addpath(stubs); rehash
cleanup = onCleanup(@() restore(stubs));

if nargin < 1 || isempty(ids), ids = 1:numel(S); end
results = struct('name',{},'status',{},'seconds',{},'message',{});
for k = ids
    fprintf('\n==================== C%02d %s ====================\n', k, S(k).name);
    tic; status = 'PASS'; msg = '';
    try
        runsection(S(k).code, dpath);
    catch ME
        status = 'FAIL';
        msg = sprintf('%s: %s', ME.identifier, ME.message);
        for s = 1:min(5,numel(ME.stack))
            msg = sprintf('%s\n    at %s:%d', msg, ME.stack(s).name, ME.stack(s).line);
        end
    end
    results(end+1) = struct('name',S(k).name,'status',status,'seconds',toc,'message',msg); %#ok<AGROW>
    close all force
end

fprintf('\n==================== SUMMARY ====================\n');
for k = 1:numel(results)
    fprintf('%-60s %s  %5.0f s\n', results(k).name, results(k).status, results(k).seconds);
    if ~isempty(results(k).message), fprintf('    %s\n', strrep(results(k).message, newline, [newline '    '])); end
end
fprintf('%d/%d passed\n', sum(strcmp({results.status},'PASS')), numel(results));
end

% -------------------------------------------------------------------------
function runsection(code, dpath) %#ok<INUSD>
% Run one section in its own workspace (dpath is used by the code)
EEG = [];
eval(code);
assert(isfield(EEG,'brainbeats'), 'No BrainBeats output')
end

function restore(stubs)
rmpath(stubs);
if isappdata(0,'brainbeats_real'), rmappdata(0,'brainbeats_real'); end
end

function writestub(file, name)
% Write a dialog function that prints its message instead of opening a window
fid = fopen(file,'w');
fprintf(fid, 'function h = %s(msg, varargin)\n', name);
fprintf(fid, '%% Test stub: print instead of opening a dialog\n');
fprintf(fid, 'fprintf(''[%s] %%s\\n'', strtrim(regexprep(char(string(msg)),''\\s+'','' '')));\n', name);
fprintf(fid, 'h = [];\n');
fclose(fid);
end

function writewrapper(file)
% Write a brainbeats_process that calls the real one with plots and gong off
L = {
    'function varargout = brainbeats_process(EEG, varargin)'
    '% Test wrapper: the repository''s brainbeats_process, plots and gong off'
    'f = getappdata(0,''brainbeats_real'');'
    'for k = {''vis_cleaning'' ''vis_outputs'' ''gong''}'
    '    i = find(strcmpi(varargin(1:2:end),k{1}))*2 - 1;'
    '    if isempty(i), varargin(end+1:end+2) = {k{1}, 0}; else, varargin{i(1)+1} = 0; end'
    'end'
    '[varargout{1:max(nargout,1)}] = f(EEG, varargin{:});'
    };
fid = fopen(file,'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
end
