% TEST_CLEAN_RR - Effect of clean_rr on beats already detected by the Muse.
%
% Loads the Muse 2 beat detections in tests/data (test_beat_detections.csv,
% with test_eeg.csv for the recording length), converts them to RR
% intervals at 250 Hz, runs clean_rr on them (no amplitudes, so the
% amplitude check is off) and plots the IBI time series and distribution
% before/after. Beats detected elsewhere are already clean, which is why
% 'rr' mode ('beat_latencies') skips clean_rr in brainbeats_process.
%
% Usage: run the script from MATLAB with EEGLAB installed. Figures are
% saved in tests/figures (clean_rr_ibi_timeseries.png,
% clean_rr_ibi_distribution.png).

clear; close all;

%% Setup (EEGLAB must be on the MATLAB path)
if ~exist('eeglab', 'file')
    error('EEGLAB not found: add its folder to the MATLAB path (addpath(''path/to/eeglab'')) and run again.')
end
eeglab nogui;

% Add BrainBeats functions to path (clean_rr lives there)
bbPath = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(bbPath, 'functions'));

testDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(testDir, 'data');
figDir  = fullfile(testDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% Load beat detections (beat times in s from the optics start)
beats = readtable(fullfile(dataDir, 'test_beat_detections.csv'));
beat_idx = beats.hrBeat == 1;
beat_sec = beats.timestamp_sec(beat_idx);
fprintf('Pre-detected beats: %d (mean IBI=%.0f ms)\n', length(beat_sec), mean(diff(beat_sec))*1000);

%% Convert to RR intervals
target_srate = 250;

% Load EEG to get the recording length in samples (sorted, unique timestamps)
eeg_raw = readtable(fullfile(dataDir, 'test_eeg.csv'));
[eeg_ts, ~] = sort(eeg_raw.ts);
[eeg_ts, ~] = unique(eeg_ts);
t_start = eeg_ts(1);
t_end   = eeg_ts(end);
npnts = length(t_start:(1/target_srate):t_end);

% Beat times to samples, keeping beats inside the recording
Rpeaks = round(beat_sec * target_srate);
Rpeaks(Rpeaks < 1) = [];
Rpeaks(Rpeaks > npnts) = [];

Rpeaks = Rpeaks(:)';  % ensure row vector
RR   = diff(Rpeaks) / target_srate;       % seconds
RR_t = Rpeaks(2:end) / target_srate;      % seconds (beat ending each interval, as get_RR)
fprintf('RR intervals: mean=%.0f ms, std=%.0f ms, n=%d\n', mean(RR)*1000, std(RR)*1000, length(RR));

%% Run clean_rr (no amplitudes: amplitude check off)
[NN, NN_t, nPeaks, idx_bad, idx_interp] = clean_rr(RR_t, RR, [], Rpeaks(2:end));
nRem = length(RR) - (length(NN) - sum(idx_interp));
idx_rem = true(1, nRem);   % for the summary count (removed beats are absent from the output)
fprintf('clean_rr: removed %d, interpolated %d, unfilled gaps %d out of %d intervals\n', ...
    nRem, sum(idx_interp), sum(idx_bad), length(RR));

% Peaks after clean_rr, without the inserted beats (NaN when no signal is given)
Rpeaks_after = [Rpeaks(1) nPeaks(~idx_interp)'];

%% IBI series
ibi_raw   = diff(beat_sec) * 1000;          % ms
ibi_after = diff(Rpeaks_after) / target_srate * 1000;  % ms

%% Figure 1: IBI time series
fig1 = figure('Position', [100 100 1200 400], 'Color', 'w');
plot(beat_sec(2:end), ibi_raw, 'g.-', 'MarkerSize', 4); hold on;
plot(Rpeaks_after(2:end)/target_srate, ibi_after, 'r.-', 'MarkerSize', 4);
xlabel('Time (s)'); ylabel('IBI (ms)');
title('clean\_rr distortion: IBI time series');
legend({'Raw beat\_detections (input)', 'After clean\_rr (output)'}, 'Location', 'best');
ylim([400 1200]);
set(gca, 'FontSize', 11);
saveas(fig1, fullfile(figDir, 'clean_rr_ibi_timeseries.png'));
fprintf('Saved: clean_rr_ibi_timeseries.png\n');

%% Figure 2: IBI distribution
fig2 = figure('Position', [100 100 600 400], 'Color', 'w');
histogram(ibi_raw, 40, 'FaceColor', [0.6 0.9 0.6], 'FaceAlpha', 0.5); hold on;
histogram(ibi_after, 40, 'FaceColor', [0.9 0.4 0.4], 'FaceAlpha', 0.5);
xlabel('IBI (ms)'); ylabel('Count');
title('clean\_rr distortion: IBI distribution');
legend({'Raw beat\_detections (input)', 'After clean\_rr (output)'}, 'Location', 'best');
set(gca, 'FontSize', 11);
saveas(fig2, fullfile(figDir, 'clean_rr_ibi_distribution.png'));
fprintf('Saved: clean_rr_ibi_distribution.png\n');

%% Summary
fprintf('\n=== SUMMARY ===\n');
fprintf('Input:  %d beats, %d RR intervals, mean=%.0f ms, std=%.0f ms\n', ...
    length(beat_sec), length(RR), mean(RR)*1000, std(RR)*1000);
fprintf('Output: %d peaks (was %d), %d IBIs\n', ...
    length(Rpeaks_after), length(Rpeaks), length(ibi_after));
fprintf('clean_rr removed %d + interpolated %d = %d flagged out of %d\n', ...
    sum(idx_rem), sum(idx_interp), sum(idx_rem)+sum(idx_interp), length(RR));
fprintf('Test complete.\n');
