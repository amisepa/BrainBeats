function results = run_tutorial_headless(ids)
% RUN_TUTORIAL_HEADLESS - Headless regression run of brainbeats_tutorial.
%
% Runs every brainbeats_tutorial call, plus 'rr' mode (HEP and HRV from beat
% latencies), a 2-ECG-channel HEP and a features run with EEG features off,
% without plots, and prints PASS/FAIL with the key outputs of each test.
% Run before a release.
%
% Usage:
%   results = run_tutorial_headless;        % all tests (EEGLAB on the path)
%   results = run_tutorial_headless(1:5);   % a subset (indices in tests())
%
% Input:
%   ids     - indices of the tests to run (default: all)
%
% Output:
%   results - struct array with fields name, status ('PASS'/'FAIL'/'SKIP'),
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
if ~exist('pop_loadset','file'), eeglab; close; end
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
% High-passed copy (EEG channels, 0.5 Hz) for the HEP runs without EEG
% cleaning: unfiltered EEG is refused there
pop_saveset(highpass_eeg(pop_loadset('filename','dataset.set','filepath',dpath)), ...
    'filename','dataset_hp.set','filepath',dpath);
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
        if isempty(EEG)
            status = 'SKIP';
        else
            assert(isfield(EEG,'brainbeats'), 'No BrainBeats output')
            info = summarize(EEG);
        end
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
fprintf('%d/%d passed, %d skipped\n', sum(strcmp({results.status},'PASS')), numel(results), sum(strcmp({results.status},'SKIP')));
end

% -------------------------------------------------------------------------
function T = tests(dpath)
% List of tests: name and a function handle returning the output EEG.
% L() reloads the scratch copy of the sample dataset, H() its high-passed
% copy (HEP without EEG cleaning).
L = @() pop_loadset('filename','dataset.set','filepath',dpath);
H = @() pop_loadset('filename','dataset_hp.set','filepath',dpath);
T = struct('name',{},'fn',{});
% brainbeats_tutorial calls, in tutorial order (plots off); its HEP on ICs
% call is 'hep_ics' below
T(end+1) = struct('name','hep_ecg','fn', @() hep_ecg_main(L));
T(end+1) = struct('name','hep_ppg','fn', @() bbp(L(),'analysis','hep','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'clean_eeg',true,'linenoise',50, ...
    'ref','infinity','highpass',1,'lowpass',20,'filttype','causal', ...
    'detectMethod','median','icamethod',1,'save',false));
T(end+1) = struct('name','feat_ecg_clean','fn', @() bbp(L(),'analysis','features','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true,'linenoise',50,'parpool',true));
T(end+1) = struct('name','feat_ppg_noclean','fn', @() bbp(pop_select(L(),'nochannel',{'ECG'}), ...
    'analysis','features','heart_signal','PPG','heart_channels',{'PPG'},'clean_eeg',false, ...
    'hrv_features', {'time' 'frequency' 'nonlinear'},'hrv_spec','welch', ...
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
T(end+1) = struct('name','coh_ecg','fn', @() bbp(pop_select(L(),'nochannel',{'PPG'}), ...
    'analysis','coherence','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',0,'parpool',0));
T(end+1) = struct('name','coh_ppg','fn', @() bbp(pop_select(L(),'nochannel',{'ECG'}), ...
    'analysis','coherence','heart_signal','PPG', ...
    'heart_channels',{'PPG'},'clean_eeg',true,'linenoise',50, ...
    'ref','infinity','ica_method',1,'parpool',false));
% Extras not in the tutorial
T(end+1) = struct('name','hep_rr_mode','fn', @() hep_rr(L, H));
T(end+1) = struct('name','rr_features','fn', @() rr_feat(L));
T(end+1) = struct('name','hep_ecg_2chan','fn', @() hep_2chan(H));
T(end+1) = struct('name','hrv_eegfeatures_off','fn', @() bbp(pop_select(L(),'nochannel',{'PPG'}), ...
    'analysis','features','heart_signal','ECG','heart_channels',{'ECG'},'eeg_features','off','parpool','off'));
T(end+1) = struct('name','hep_baseline_export','fn', @() hep_blreg(H));
T(end+1) = struct('name','hep_adaptive_window','fn', @() bbp(pop_select(H(),'nochannel',{'PPG'}), ...
    'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'},'hep_window','adaptive','clean_eeg',false,'save',0));
T(end+1) = struct('name','hep_ppg_transit','fn', @() hep_ppg_transit(H));
T(end+1) = struct('name','hep_tf_surrogates','fn', @() hep_tf_surr(H));
T(end+1) = struct('name','hep_ics','fn', @() hep_ics(L, true));
T(end+1) = struct('name','hep_ics_noclean','fn', @() hep_ics(H, false));
T(end+1) = struct('name','hep_csd','fn', @() bbp(pop_select(L(),'nochannel',{'PPG'}),'analysis','hep', ...
    'heart_signal','ECG','heart_channels',{'ECG'},'clean_eeg',true,'icamethod',1,'ref','csd','save',0));
T(end+1) = struct('name','hep_ppg_pat_eeg','fn', @() hep_ppg_pat_eeg(H));
T(end+1) = struct('name','feat_csd','fn', @() bbp(pop_select(L(),'nochannel',{'PPG'}),'analysis','features', ...
    'heart_signal','ECG','heart_channels',{'ECG'},'clean_eeg',true,'icamethod',1,'ref','csd', ...
    'eeg_features',{'frequency'},'hrv_features','off','save',0));
T(end+1) = struct('name','hep_gedai','fn', @() hep_gedai(L));
T(end+1) = struct('name','hep_ecg_regression','fn', @() hep_ecgreg(L));
T(end+1) = struct('name','input_errors','fn', @() input_errors(L, H));
T(end+1) = struct('name','feat_relative_power','fn', @() feat_relative(L));
T(end+1) = struct('name','hep_ica_ecg','fn', @() hep_removal(L, 'ica_ecg'));
T(end+1) = struct('name','hep_ssp','fn', @() hep_removal(L, 'ssp'));
T(end+1) = struct('name','hep_ica_caracas','fn', @() hep_removal(L, 'ica_caracas'));
T(end+1) = struct('name','hep_ica_epochs_60hz','fn', @() hep_ica_opts(L, 'epochs'));
T(end+1) = struct('name','hep_ica_lowpass_noasr','fn', @() hep_ica_opts(L, 'lowpass'));
end

function EEG = bbp(EEG, varargin)
% brainbeats_process with plots and gong forced off
for k = {'vis_cleaning' 'vis_outputs' 'gong'}
    i = find(strcmpi(varargin,k{1}));
    if isempty(i), varargin(end+1:end+2) = {k{1}, 0}; else, varargin{i+1} = 0; end
end
EEG = brainbeats_process(EEG, varargin{:});
end

function EEG = hep_ecg_main(L)
% The brainbeats_tutorial call: ICA fitted on the continuous copy, whose line
% noise (50 Hz in the sample data, 'linenoise' not given) was detected and
% notched
EEG = bbp(L(),'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'},'clean_eeg',true, ...
    'ica_method',1,'hep_surrogates',100,'keep_heart',true);
I = EEG.brainbeats.preprocessings.ica;
assert(~I.epoched && I.copy_options && isequal(I.line_noise_notched, 50), ...
    'ICA copy: epoched %d, copy options %d, line noise notched %s', I.epoched, I.copy_options, mat2str(I.line_noise_notched))
fprintf('ICA on the continuous copy (%d samples, rank %d), line noise notched at %g Hz\n', I.samples, I.rank, I.line_noise_notched);
end

function EEG = hep_rr(L, H)
% Beats from an HRV-only run, then HEP from those latencies ('rr' mode)
E = bbp(pop_select(L(),'channel',{'ECG'}),'analysis','features','eeg','off', ...
    'heart_signal','ECG','heart_channels',{'ECG'},'save',0);
EEG = bbp(pop_select(H(),'nochannel',{'ECG' 'PPG'}),'analysis','hep','heart_signal','rr', ...
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
% HRSP/HRPC of all channels and 50 surrogates: the HEP rebuilt from the
% continuous data must equal the epochs, and the stats must be complete
EEG = bbp(pop_select(L(),'nochannel',{'PPG'}),'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'hep_surrogates',50,'save',0);
tf = EEG.brainbeats.hrsp; S = EEG.brainbeats.surrogate;
iE = ismember({EEG.chanlocs.labels}, S.channels);
assert(max(abs(mean(EEG.data(iE,:,:),3) - S.hep.real),[],'all') < 1e-6, 'HEP from the continuous data differs from the epochs')
assert(isequal(size(tf.hrsp), size(tf.hrpc), size(S.hrsp.p_fdr), size(tf.hrsp_induced), size(S.hrsp_induced.p_fdr)) ...
    && ~any(isnan(tf.hrsp(:))), 'Incomplete HRSP/HRPC outputs')
fprintf('HRSP/HRPC: %d channels x %d freqs x %d times; %.1f%% of HEP points differ from the surrogates (FDR)\n', ...
    size(tf.hrsp), 100*mean(S.hep.p_fdr(:) < .05));
% Pre-R-peak values of each heartbeat and the QRS map
P = EEG.brainbeats.prestim; C = EEG.brainbeats.cfa; nE = sum(iE);
assert(isequal(size(P.baseline), [nE EEG.trials]) && isequal(size(P.power), [nE 3 EEG.trials]) && ...
    mean(isnan(P.power(:))) < .05 && numel(P.beat_time) == EEG.trials, 'Incomplete pre-R-peak values')
assert(numel(C.map) == nE && abs(C.latency) <= 60 && ~isempty(C.notes), 'Incomplete QRS map')
fprintf('Pre-R-peak values: %d channels x %d bands x %d heartbeats; QRS map at %g ms, %s\n', size(P.power), C.latency, C.summary);
end

function EEG = hep_blreg(L)
% Exported baselines (no correction of the HEP): each heartbeat's mean in
% 'hep_baseline_win'; the removed 'hep_baseline' option only warns; the
% regression with BASELINE_REGRESSION keeps the average and can be undone
args = {'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'}, ...
    'clean_eeg',false,'keep_heart',true,'hep_surrogates',0,'save',0};
EEG = bbp(pop_select(L(),'nochannel',{'PPG'}), args{:}, 'hep_baseline_win',[-150 -50]);
P = EEG.brainbeats.prestim;
iE = ismember({EEG.chanlocs.labels}, P.channels);
iB = EEG.times >= -150 & EEG.times <= -50;
assert(max(abs(P.baseline - squeeze(mean(double(EEG.data(iE,iB,:)),2))),[],'all') < 1e-6, 'Exported baselines differ from the epochs')
% 'keep_heart': the heart channel is filtered (raw sample ECG: epoch offset 1.1
% times the range of its average; filtered: 0.007)
h = squeeze(double(EEG.data(strcmp({EEG.chanlocs.labels},'ECG'),:,:)));
r = median(abs(mean(h,1))) / range(mean(h,2));
assert(r < 0.1, 'Kept heart channel not filtered: epoch offset %.2f times the range of its average', r)
s = warning('error', 'BrainBeats:removed');
try
    bbp(pop_select(L(),'nochannel',{'PPG'}), args{:}, 'hep_baseline','regression');
    ok = false;
catch ME
    ok = strcmp(ME.identifier, 'BrainBeats:removed');
end
warning(s);
assert(ok, '''hep_baseline'',''regression'' did not warn that it was removed')
[Y, beta, bl] = baseline_regression(EEG.data(iE,:,:), EEG.times, [-150 -50]);
assert(max(abs(mean(Y,3) - mean(EEG.data(iE,:,:),3)),[],'all') < 1e-6, 'Average HEP changed')
undo = Y + beta .* permute(bl - mean(bl,2), [1 3 2]);
assert(max(abs(undo - EEG.data(iE,:,:)),[],'all') < 1e-6, 'Correction cannot be undone')
end

function EEG = hep_ics(L, clean)
% HEP measures of all ICs (ICA of the cleaning: the brainbeats_tutorial call;
% or run for them when the EEG is not cleaned): one row per IC, with maps and
% ICLabel classes
if clean
    EEG = bbp(pop_select(L(),'nochannel',{'PPG'}),'analysis','hep','heart_signal','ECG', ...
        'heart_channels',{'ECG'},'clean_eeg',true,'ica_method',1, ...
        'hep_level','both','hep_surrogates',100);
else
    EEG = bbp(pop_select(L(),'nochannel',{'PPG'}),'analysis','hep','heart_signal','ECG', ...
        'heart_channels',{'ECG'},'clean_eeg',false,'icamethod',1,'hep_level','both','hep_surrogates',20,'save',0);
end
I = EEG.brainbeats.ics;
nIC = numel(I.ic);
assert(nIC > 1 && size(I.hep,1) == nIC && size(I.tf.hrsp,1) == nIC && size(I.maps,2) == nIC, 'Incomplete IC outputs')
assert(isequal(size(I.tf.hrsp), size(I.surrogate.hrsp.p_fdr)) && size(I.iclabel,1) == nIC, 'Incomplete IC statistics')
assert(isfield(EEG.brainbeats,'roi'), 'ROI outputs missing')
assert(isempty(intersect(I.ic, I.removed_ics)) && max([I.ic; I.removed_ics]) == numel(I.ic) + numel(I.removed_ics), 'IC numbers do not match the decomposition')
fprintf('%d ICs: HEP, HRSP, HRPC and surrogates; %d HEP points FDR-significant\n', nIC, sum(I.surrogate.hep.p_fdr(:) < .05));
end

function EEG = hep_ppg_pat_eeg(L)
% PPG without ECG: no pulse arrival time unless one is given ('auto' is an
% error); with 'eeg' it comes from the EEG cardiac field artifact and must be
% close to the ECG-based one (112 ms to the pulse onsets on the sample)
args = {'analysis','hep','heart_signal','PPG','heart_channels',{'PPG'},'clean_eeg',false,'save',0};
id = '';
try
    bbp(pop_select(L(),'nochannel',{'ECG'}), args{:});
catch ME
    id = ME.identifier;
end
assert(strcmp(id,'BrainBeats:input'), 'PPG without ECG and without a pulse arrival time did not raise an error')
EEG = bbp(pop_select(L(),'nochannel',{'ECG'}), args{:}, 'ppg_transit','eeg');
p = EEG.brainbeats.preprocessings.ppg_transit;
assert(contains(p.method,'EEG') && abs(p.pat - 112) < 50, 'PAT from the EEG: %g ms (%s)', p.pat, p.method)
fprintf('PAT from the EEG: %.0f ms (ECG: 112 ms)\n', p.pat);
end

function EEG = hep_ecgreg(L)
% Heart artifact removed by ECG regression: the heart-locked EEG amplitude
% must drop by at least half, with little of the EEG variance removed, and
% no heart component removed by ICA
EEG = bbp(pop_select(L(),'nochannel',{'PPG'}),'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true,'icamethod',1,'heart_removal','ecg_regression','save',0);
r = EEG.brainbeats.preprocessings.heart_regression;
assert(r.cfa_after < 0.5*r.cfa_before && r.variance_removed < 5, ...
    'ECG regression: %.2f -> %.2f uV, %.1f%% of the variance removed', r.cfa_before, r.cfa_after, r.variance_removed)
fprintf('ECG regression: %.2f -> %.2f uV (%.0f%% reduction), %.1f%% of the variance removed\n', ...
    r.cfa_before, r.cfa_after, 100*(1 - r.cfa_after/r.cfa_before), r.variance_removed);
end

function EEG = hep_removal(L, method)
% Heart artifact removal options: they run and report the heart-locked
% amplitude before and after
EEG = bbp(pop_select(L(),'nochannel',{'PPG'}),'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true,'icamethod',1,'heart_removal',method,'save',0);
pre = EEG.brainbeats.preprocessings;
if strcmp(method,'ssp'), r = pre.heart_ssp; else, r = pre; end
assert(isfinite(r.cfa_before) && isfinite(r.cfa_after) && isfinite(r.cfa_late_after), 'Heart-locked amplitude missing')
if strcmp(method,'ssp')
    % the projected topographies stay out after the ICA cleaning (the ICA
    % copy is projected too), on the channels the projection was computed
    % on (the bad channels are interpolated afterwards)
    lab = {EEG.chanlocs.labels};
    iE = find(~ismember(lab, [pre.removed_eeg_channels {'ECG' 'PPG'}]));
    assert(numel(iE) == size(r.patterns,1), 'SSP: channel count differs from the projection')
    Y = reshape(double(EEG.data(iE,:,:)), numel(iE), []);
    left = norm(r.patterns' * Y, 'fro')^2 / norm(Y, 'fro')^2;
    assert(left < 1e-6, 'SSP: %.2g%% of the output EEG energy along the projected topographies', 100*left)
end
fprintf('%s: %.2f -> %.2f uV (-50 to 100 ms)\n', method, r.cfa_before, r.cfa_after);
end

function EEG = hep_ica_opts(L, which)
% Options of the copy the HEP ICA is fitted on.
% 'epochs': 60-Hz line noise (the sample's 50-Hz line noise shifted to 60 Hz,
% so that it keeps its distribution over channels: a sinusoid common to all
% channels is removed by the average reference), the copy epoched and
% resampled to 125 Hz, and Picard made to diverge (a stub returning NaN
% weights): the copy must be notched at 60 Hz, epoched, at 125 Hz, and the
% ICA rerun with Infomax.
% 'lowpass': copy low-passed at 40 Hz (no line-noise notch then) without ASR.
EEG = pop_select(L(),'nochannel',{'PPG'});
args = {'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'},'clean_eeg',true,'icamethod',1,'save',0};
if strcmp(which, 'epochs')
    iE = ~strcmp({EEG.chanlocs.labels}, 'ECG');
    F = pop_eegfiltnew(EEG,'locutoff',48,'hicutoff',52,'revfilt',1);     % 50 Hz notched out
    B = pop_eegfiltnew(EEG,'locutoff',49.5,'hicutoff',50.5,'filtorder',8*round(EEG.srate));   % the 50-Hz line only
    t = (0:EEG.pnts-1) / EEG.srate;
    shifted = real(hilbert(double(B.data(iE,:))') .* exp(1i*2*pi*10*t'))';   % moved to 60 Hz
    EEG.data(iE,:) = F.data(iE,:) + shifted;
    stub = fullfile(tempdir, 'brainbeats_test_picard');
    if ~exist(stub,'dir'), mkdir(stub); end
    fid = fopen(fullfile(stub,'picard.m'), 'w');
    fprintf(fid, 'function [Y, W] = picard(X, varargin)\n%% Test stub: a diverged Picard (NaN weights)\nY = X; W = nan(size(X,1));\n');
    fclose(fid);
    addpath(stub); c = onCleanup(@() rmpath(stub));
    EEG = bbp(EEG, args{:}, 'ica_epochs',true, 'ica_srate',125);
    I = EEG.brainbeats.preprocessings.ica;
    assert(isequal(I.line_noise_notched, 60) && I.epoched && I.srate == 125 && I.method == 2, ...
        'ICA copy: line noise %s Hz, epoched %d, %g Hz, method %d (expected 60 Hz, epoched, 125 Hz, Infomax after Picard failed)', ...
        mat2str(I.line_noise_notched), I.epoched, I.srate, I.method)
    assert(all(isfinite(EEG.icaweights(:))), 'Non-finite ICA weights after the fallback')
else
    EEG = bbp(EEG, args{:}, 'ica_lowpass',40, 'ica_asr',0);
    I = EEG.brainbeats.preprocessings.ica;
    P = EEG.brainbeats.parameters;
    assert(I.copy_options && ~I.epoched && isempty(I.line_noise_notched) && P.ica_asr == 0 && P.ica_lowpass == 40, ...
        'ICA copy: options %d, epoched %d, line noise %s, ica_asr %g, ica_lowpass %g', I.copy_options, I.epoched, ...
        mat2str(I.line_noise_notched), P.ica_asr, P.ica_lowpass)
end
fprintf('ICA copy (%s): %d samples at %g Hz, epoched %d, line noise notched %s, method %d, %d ICs removed\n', which, ...
    I.samples, I.srate, I.epoched, mat2str(I.line_noise_notched), I.method, numel(EEG.brainbeats.preprocessings.removed_eeg_components));
end

function EEG = hep_gedai(L)
% GEDAI cleaning, with the IC measures (only if the plugin is installed: not
% installed by the test)
if ~exist('GEDAI','file')
    fprintf('SKIPPED: GEDAI plugin not installed.\n');
    EEG = []; return
end
EEG = bbp(pop_select(L(),'nochannel',{'PPG'}),'analysis','hep','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',true,'clean_method','gedai','hep_level','ics','icamethod',1,'save',0);
assert(isfield(EEG.brainbeats,'ics') && isfield(EEG.etc,'GEDAI'), 'GEDAI or IC outputs missing')
assert(isfield(EEG.brainbeats.preprocessings,'cfa_gedai'), 'Heart-locked amplitude of GEDAI missing')
end

function EEG = input_errors(L, H)
% On the command line, failed checks and unknown inputs must raise errors
EEG = pop_select(H(),'nochannel',{'PPG'});
good = {'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'},'clean_eeg',false,'save',0};
% Unfiltered EEG for HEP without cleaning (the HEP is not baseline-corrected),
% and sampling rates above 1024 Hz, are refused; above 512 Hz, a warning
E = pop_select(L(),'nochannel',{'PPG'});
E2 = EEG; E2.srate = 2048;
E3 = EEG; E3.srate = 1000;
s = warning('error', 'BrainBeats:srate');
cases = {E, 'BrainBeats:input', 'unfiltered EEG'; E2, 'BrainBeats:input', '2048 Hz'; E3, 'BrainBeats:srate', '1000 Hz'};
for k = 1:size(cases,1)
    id = ''; txt = '';
    try
        bbp(cases{k,1}, good{:});
    catch ME
        id = ME.identifier; txt = ME.message;
    end
    assert(strcmp(id, cases{k,2}), '%s: expected %s, got ''%s''', cases{k,3}, cases{k,2}, id)
    fprintf('%s refused: %s\n', cases{k,3}, txt(1:min(end,90)));
end
warning(s);
bad = { [good(1:4) {'heart_channels',{'NOPE'}} good(7:end)], ...   % heart channel not in the data
    [good {'lowpas',20}], ...                                      % typo in an input
    [good {'vis_outputs'}] };                                      % not 'key', value pairs
for k = 1:numel(bad)
    id = '';
    try
        bbp(EEG, bad{k}{:});
    catch ME
        id = ME.identifier;
    end
    assert(strcmp(id,'BrainBeats:input'), 'Bad input %d did not raise a BrainBeats error (got ''%s'')', k, id)
end
% A value equal to a key is not read as that key
P = getparams_cmd('analysis','hep','heart_signal','ppg','heart_channels',{'PPG'},'ppg_transit','eeg', ...
    'clean_eeg',false,'parpool',0,'gpu',0);
assert(isequal(P.eeg,true) && strcmp(P.ppg_transit,'eeg'), '''ppg_transit'',''eeg'' was read as the ''eeg'' input')
% Switches accept 'on'/'off' ('clean_eeg','off' cleaned and 'save','off' saved before v1.7)
P = getparams_cmd('analysis','hep','heart_signal','ecg','heart_channels',{'ECG'},'clean_eeg','off', ...
    'save','off','hep_cfa','off','vis_outputs','on','keep_heart','off','parpool',0,'gpu',0);
assert(isequal([P.clean_eeg P.save P.hep_cfa P.vis_outputs P.keep_heart], [false false false true false]), '''on''/''off'' switches misread')
% Out-of-range values raise errors
bad = {{'clean_eeg','maybe'}, {'clean_eeg',true,'hep_baseline_win',[100 200]}, ...
    {'clean_eeg',true,'ica_lowpass',0.5}, {'clean_eeg',true,'ica_srate',5}};
for k = 1:numel(bad)
    failed = false;
    try
        getparams_cmd('analysis','hep','heart_signal','ecg','heart_channels',{'ECG'},bad{k}{:},'parpool',0,'gpu',0);
    catch
        failed = true;
    end
    assert(failed, 'No error for %s', strjoin(cellfun(@(v) char(string(mat2str(v))), bad{k}, 'UniformOutput', false), ' '))
end
EEG = bbp(EEG, good{:});
% A dataset without a file is saved under its name
[name, folder] = output_file(struct('filename','','filepath','','setname','my data 1'), '_HEP.set');
assert(strcmp(name,'my_data_1_HEP.set') && strcmp(folder,pwd), 'Output file: %s', name)
end

function EEG = feat_relative(L)
% 'eeg_norm',2: the band values are the relative powers, unrounded
EEG = bbp(pop_select(L(),'nochannel',{'PPG'}),'analysis','features','eeg','on','heart_signal','ECG', ...
    'heart_channels',{'ECG'},'clean_eeg',false,'hrv_features','off','eeg_features',{'frequency'},'eeg_norm',2,'save',0);
F = EEG.brainbeats.features.EEG;
for b = {'delta' 'theta' 'alpha' 'beta' 'gamma'}
    assert(isequal(F.frequency.(b{1}), F.qeeg.([b{1} '_rel'])), '%s is not the relative power', b{1})
end
tot = F.frequency.delta + F.frequency.theta + F.frequency.alpha + F.frequency.beta + F.frequency.gamma;
assert(all(tot > 0.5 & tot <= 1.01), 'Relative powers sum to %.2f-%.2f', min(tot), max(tot))
assert(any(mod(F.frequency.gamma*1e3, 1) ~= 0), 'Band values are rounded')
fprintf('Relative power: bands sum to %.2f-%.2f of the total; gamma %.5f on average\n', min(tot), max(tot), mean(F.frequency.gamma));
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
