function results = run_tutorial_headless(ids)
% RUN_TUTORIAL_HEADLESS - Headless regression run of brainbeats_tutorial.
%
% Runs every brainbeats_tutorial call, plus 'rr' mode (HEP and HRV from beat
% latencies), a 2-ECG-channel HEP and a features run with EEG features off,
% without plots, and prints PASS/FAIL with the key outputs of each test.
% Run before a release.
%
% Usage:
%   matlab -batch "addpath('path/to/eeglab'); cd('path/to/BrainBeats/tests'); run_tutorial_headless"
%   results = run_tutorial_headless;        % all tests
%   results = run_tutorial_headless(1:5);   % a subset (indices in tests())
%
% Input:
%   ids     - indices of the tests to run (default: all)
%
% Output:
%   results - struct array with fields name, status ('PASS'/'FAIL'),
%             seconds, message (error and stack) and summary (key outputs)
%
% Notes:
% - Plots and the end gong are forced off ('vis_cleaning', 'vis_outputs' and
%   'gong' = 0): drawing can hang under 'matlab -batch', so the plotting
%   code still needs one interactive run of brainbeats_tutorial.
% - errordlg/warndlg are replaced by stubs that print the message, and the
%   sample data are copied to tempdir so 'save' never writes into the repo.
% - Any other copy of BrainBeats on the path (e.g. eeglab/plugins/BrainBeats*)
%   is removed so this repository's code is what runs.
% - Tests with 'clean_eeg' and 'ref','infinity' need the REST plugin, which
%   BrainBeats installs if missing.
%
% Runtime: ~10-20 min (ICA and nonlinear EEG features dominate).
%
% Copyright (C) - Cedric Cannard, 2026

% Paths: EEGLAB, then this repository instead of any other BrainBeats copy
bbroot = fileparts(fileparts(mfilename('fullpath')));
if ~exist('pop_loadset','file'), eeglab nogui; end
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bbroot));
if ~isempty(old), rmpath(old{:}); end
addpath(bbroot, fullfile(bbroot,'functions'));

% Dialog stubs (printed instead of opened)
stubs = fullfile(tempdir,'brainbeats_test_stubs');
if ~exist(stubs,'dir'), mkdir(stubs); end
writestub(fullfile(stubs,'errordlg.m'), 'errordlg');
writestub(fullfile(stubs,'warndlg.m'),  'warndlg');
addpath(stubs);
cleanup = onCleanup(@() rmpath(stubs));

% Scratch copy of the sample data
dpath = fullfile(tempdir,'brainbeats_test_data');
if ~exist(dpath,'dir'), mkdir(dpath); end
copyfile(fullfile(bbroot,'sample_data','dataset.set'), dpath);
fprintf('Testing %s\n', which('brainbeats_process'));

% Run the tests; errors are caught and reported, not thrown
T = tests(dpath);
if nargin < 1 || isempty(ids), ids = 1:numel(T); end
results = struct('name',{},'status',{},'seconds',{},'message',{},'summary',{});
for k = ids
    fprintf('\n==================== T%02d %s ====================\n', k, T(k).name);
    tic; status = 'PASS'; msg = ''; info = '';
    try
        EEG = T(k).fn();
        info = summarize(EEG);
    catch ME
        status = 'FAIL';
        msg = sprintf('%s: %s', ME.identifier, ME.message);
        for s = 1:min(5,numel(ME.stack))
            msg = sprintf('%s\n    at %s:%d', msg, ME.stack(s).name, ME.stack(s).line);
        end
    end
    results(end+1) = struct('name',T(k).name,'status',status,'seconds',toc,'message',msg,'summary',info); %#ok<AGROW>
    close all force
end

fprintf('\n==================== SUMMARY ====================\n');
for k = 1:numel(results)
    fprintf('%-22s %s  %5.0f s\n', results(k).name, results(k).status, results(k).seconds);
    if ~isempty(results(k).message), fprintf('    %s\n', strrep(results(k).message, newline, [newline '    '])); end
    fprintf('%s', results(k).summary);
end
fprintf('%d/%d passed\n', sum(strcmp({results.status},'PASS')), numel(results));
end

% -------------------------------------------------------------------------
function T = tests(dpath)
% List of tests: name and a function handle returning the output EEG.
% L() reloads the scratch copy of the sample dataset.
L = @() pop_loadset('filename','dataset.set','filepath',dpath);
T = struct('name',{},'fn',{});
% brainbeats_tutorial calls, in tutorial order (plots off)
T(end+1) = struct('name','hep_ecg','fn', @() bbp(L(),'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true,'ica_method',1,'keep_heart',true));
T(end+1) = struct('name','hep_ppg','fn', @() bbp(L(),'analysis','hep','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'clean_eeg',true,'linenoise',50, ...
    'ref','infinity','highpass',.5,'lowpass',20,'filttype','causal', ...
    'detectMethod','median','icamethod',1,'save',false));
T(end+1) = struct('name','feat_ecg_clean','fn', @() bbp(L(),'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true,'linenoise',50,'parpool',true));
T(end+1) = struct('name','feat_ppg_noclean','fn', @() bbp(pop_select(L(),'nochannel',{'ECG'}), ...
    'analysis','features','heart_signal','PPG','heart_channels',{'PPG'},'clean_eeg',false,'linenoise',50, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'},'hrv_spec','LombScargle', ...
    'eeg_features', {'time' 'frequency'},'eeg_norm',0,'parpool','on','save',false));
T(end+1) = struct('name','hrv_ecg_eegoff','fn', @() bbp(pop_select(L(),'nochannel',{'PPG'}), ...
    'analysis','features','eeg','off','heart_signal','ECG','heart_channels',{'ECG'}, ...
    'hrv_features',{'time' 'frequency' 'nonlinear'}));
T(end+1) = struct('name','hrv_ecg_noeeg','fn', @() bbp(pop_select(L(),'channel',{'ECG'}), ...
    'analysis','features','eeg','off','heart_signal','ECG','heart_channels',{'ECG'}));
T(end+1) = struct('name','hrv_ppg_noeeg','fn', @() bbp(pop_select(L(),'channel',{'PPG'}), ...
    'analysis','features','eeg','off','heart_signal','PPG','heart_channels',{'PPG'}));
T(end+1) = struct('name','eeg_only','fn', @() bbp(pop_select(L(),'nochannel',{'PPG' 'ECG'}), ...
    'analysis','features','heart_signal','off', ...
    'eeg_features',{'time' 'frequency','nonlinear'},'clean_eeg',1,'linenoise',50, ...
    'ref','infinity','ica_method',1,'parpool',1,'save',1));
T(end+1) = struct('name','rm_heart','fn', @() bbp(L(), 'analysis', 'rm_heart', ...
    'heart_signal', 'ECG', 'heart_channels', {'ECG'}, ...
    'clean_eeg', true, 'filttype', 'noncausal','linenoise', 50, ...
    'conf_thresh', .65, 'ica_method', 2, 'keep_heart',true));
T(end+1) = struct('name','coh_ecg','fn', @() bbp(pop_select(L(),'nochannel',{'PPG'}), ...
    'analysis','coherence','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',0,'linenoise',50,'ref','infinity','ica_method',1,'parpool',0));
T(end+1) = struct('name','coh_ppg','fn', @() bbp(pop_select(L(),'nochannel',{'ECG'}), ...
    'analysis','coherence','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'clean_eeg',true,'linenoise',50, ...
    'ref','infinity','ica_method',1,'parpool',false));
% Extras not in the tutorial
T(end+1) = struct('name','hep_rr_mode','fn', @() hep_rr(L));
T(end+1) = struct('name','rr_features','fn', @() rr_feat(L));
T(end+1) = struct('name','hep_ecg_2chan','fn', @() hep_2chan(L));
T(end+1) = struct('name','hrv_eegfeatures_off','fn', @() bbp(pop_select(L(),'nochannel',{'PPG'}), ...
    'analysis','features','heart_signal','ECG','heart_channels',{'ECG'},'eeg_features','off','parpool','off'));
T(end+1) = struct('name','hep_baseline_regression','fn', @() hep_blreg(L));
T(end+1) = struct('name','hep_adaptive_window','fn', @() bbp(pop_select(L(),'nochannel',{'PPG'}), ...
    'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'},'hep_window','adaptive','clean_eeg',false,'save',0));
T(end+1) = struct('name','hep_ppg_transit','fn', @() hep_ppg_transit(L));
T(end+1) = struct('name','hep_tf_surrogates','fn', @() hep_tf_surr(L));
end

function EEG = bbp(EEG, varargin)
% brainbeats_process with plots and gong forced off
for k = {'vis_cleaning' 'vis_outputs' 'gong'}
    i = find(strcmpi(varargin,k{1}));
    if isempty(i), varargin(end+1:end+2) = {k{1}, 0}; else, varargin{i+1} = 0; end
end
EEG = brainbeats_process(EEG, varargin{:});
end

function EEG = hep_rr(L)
% Beats from an HRV-only run, then HEP from those latencies ('rr' mode)
E = bbp(pop_select(L(),'channel',{'ECG'}),'analysis','features','eeg','off', ...
    'heart_signal','ECG','heart_channels',{'ECG'},'save',0);
EEG = bbp(pop_select(L(),'nochannel',{'ECG' 'PPG'}),'analysis','hep','heart_signal','rr', ...
    'beat_latencies',E.brainbeats.preprocessings.NN_times,'clean_eeg',false,'save',0);
end

function EEG = rr_feat(L)
% HRV from beat latencies must match the ECG path when the beats are the same
E = bbp(pop_select(L(),'channel',{'ECG'}),'analysis','features','eeg','off', ...
    'heart_signal','ECG','heart_channels',{'ECG'},'save',0);
NN = E.brainbeats.preprocessings.NN(:);
beats = [0; cumsum(NN)] + E.brainbeats.preprocessings.NN_times(1) - NN(1);
EEG = bbp(pop_select(L(),'nochannel',{'ECG' 'PPG'}),'analysis','features','heart_signal','rr', ...
    'beat_latencies',beats,'eeg','off','save',0);
a = EEG.brainbeats.features.HRV.time.RMSSD;  b = E.brainbeats.features.HRV.time.RMSSD;
assert(abs(a-b) < 1, 'RMSSD from rr mode (%g) differs from the ECG path (%g)', a, b)
end

function EEG = hep_2chan(L)
% Two ECG channels (second = noisier copy): the cleaner one must be used
EEG = pop_select(L(),'nochannel',{'PPG'});
iE = find(strcmp({EEG.chanlocs.labels},'ECG'));
rng(1);
EEG.data(end+1,:) = EEG.data(iE,:) + 0.3*std(EEG.data(iE,:))*randn(1,EEG.pnts);
EEG.nbchan = EEG.nbchan+1;
EEG.chanlocs(end+1) = EEG.chanlocs(iE); EEG.chanlocs(end).labels = 'ECG2';
EEG = eeg_checkset(EEG);
EEG = bbp(EEG,'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG2' 'ECG'},'clean_eeg',false,'save',0);
used = EEG.brainbeats.preprocessings.heart_channel_used;
assert(strcmp(used,'ECG'), 'The noisier heart channel (%s) was selected', used)
end

function EEG = hep_ppg_transit(L)
% PPG beats shifted by the pulse arrival time estimated from the ECG: the
% HEP must then resemble the ECG-locked HEP
EEG = bbp(L(),'analysis','hep','heart_signal','PPG','heart_channels',{'PPG'}, ...
    'ppg_transit','ECG','clean_eeg',false,'save',0);
E = bbp(pop_select(L(),'nochannel',{'PPG'}),'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'save',0);
pat = EEG.brainbeats.preprocessings.ppg_transit.median;
assert(pat > 100 && pat < 600, 'Implausible pulse arrival time (%g ms)', pat)
t = E.times >= 0 & E.times < 600;
mE = mean(E.data(:,t,:),3);  mP = mean(EEG.data(:,t,:),3);
r = mean(arrayfun(@(k) corr(mE(k,:)', mP(k,:)'), 1:size(mE,1)));
assert(r > 0.7, 'PPG HEP does not match the ECG HEP after transit correction (r = %.2f)', r)
fprintf('Pulse arrival time %.0f ms; HEP correlation with the ECG HEP r = %.2f\n', pat, r);
end

function EEG = hep_tf_surr(L)
% HRSP/HEPC of all channels and 50 surrogates: the HEP rebuilt from the
% continuous data must equal the epochs, and the stats must be complete
EEG = bbp(pop_select(L(),'nochannel',{'PPG'}),'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'hep_tf',true,'hep_surrogates',50,'save',0);
tf = EEG.brainbeats.hrsp; S = EEG.brainbeats.surrogate;
iE = ismember({EEG.chanlocs.labels}, S.channels);
assert(max(abs(mean(EEG.data(iE,:,:),3) - S.hep.real),[],'all') < 1e-6, 'HEP from the continuous data differs from the epochs')
assert(isequal(size(tf.hrsp), size(tf.hepc), size(S.hrsp.p_fdr)) && ~any(isnan(tf.hrsp(:))), 'Incomplete HRSP/HEPC outputs')
fprintf('HRSP/HEPC: %d channels x %d freqs x %d times; %.1f%% of HEP points differ from the surrogates (FDR)\n', ...
    size(tf.hrsp), 100*mean(S.hep.p_fdr(:) < .05));
end

function EEG = hep_blreg(L)
% Regression baseline: same average HEP, heart channel untouched, reversible
args = {'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'}, ...
    'clean_eeg',false,'keep_heart',true,'save',0};
A = bbp(pop_select(L(),'nochannel',{'PPG'}), args{:});
EEG = bbp(pop_select(L(),'nochannel',{'PPG'}), args{:}, 'hep_baseline','regression');
r = EEG.brainbeats.preprocessings.baseline_regression;
iE = ismember({EEG.chanlocs.labels}, r.channels);
assert(max(abs(mean(A.data(iE,:,:),3) - mean(EEG.data(iE,:,:),3)),[],'all') < 1e-6, 'Average HEP changed')
assert(isequal(A.data(~iE,:,:), EEG.data(~iE,:,:)), 'Heart channel was corrected')
undo = EEG.data(iE,:,:) + r.beta .* permute(r.baseline - mean(r.baseline,2), [1 3 2]);
assert(max(abs(undo - A.data(iE,:,:)),[],'all') < 1e-6, 'Correction cannot be undone')
end

% -------------------------------------------------------------------------
function s = summarize(EEG)
% Text summary of the output: data size, preprocessing and HRV scalars,
% mean EEG band power and coherence range
s = sprintf('    data %s, trials %d, events %d\n', mat2str(size(EEG.data)), EEG.trials, numel(EEG.event));
if ~isfield(EEG,'brainbeats'), return; end
bb = EEG.brainbeats;
if isfield(bb,'preprocessings'), s = [s flat(bb.preprocessings,'pre',1)]; end
if isfield(bb,'features') && isstruct(bb.features) && isfield(bb.features,'HRV')
    s = [s flat(bb.features.HRV,'HRV',1)];
end
if isfield(bb,'features') && isstruct(bb.features) && isfield(bb.features,'EEG') && isfield(bb.features.EEG,'frequency')
    fq = bb.features.EEG.frequency;
    for b = {'delta' 'theta' 'alpha' 'beta' 'gamma'}
        if isfield(fq,b{1}), s = [s sprintf('    EEG.%s mean %.3g\n', b{1}, mean(fq.(b{1})(:),'omitnan'))]; end %#ok<AGROW>
    end
end
if isfield(bb,'coherence')
    c = bb.coherence.coh(:);
    s = [s sprintf('    coherence: range [%.3g %.3g], NaN %d\n', min(c), max(c), sum(isnan(c)))];
end
end

function s = flat(x, pre, depth)
% Print the scalar and short-string fields of a struct (up to 3 levels deep)
s = '';
if depth > 3, return; end
fn = fieldnames(x);
for i = 1:numel(fn)
    v = x.(fn{i});
    if isstruct(v) && isscalar(v)
        s = [s flat(v, [pre '.' fn{i}], depth+1)]; %#ok<AGROW>
    elseif (isnumeric(v) || islogical(v)) && isscalar(v)
        s = [s sprintf('    %s.%s = %.4g\n', pre, fn{i}, double(v))]; %#ok<AGROW>
    elseif ischar(v) && size(v,1) == 1 && numel(v) < 40
        s = [s sprintf('    %s.%s = %s\n', pre, fn{i}, v)]; %#ok<AGROW>
    end
end
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
