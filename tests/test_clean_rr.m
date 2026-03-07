% test_clean_rr.m
% Demonstrate that clean_rr distorts already-clean pre-detected beat intervals.
%
% Loads Muse 2 beat_detections, feeds them through BrainBeats clean_rr,
% and plots IBI distribution + time series before/after to show the problem.

clear; close all;

%% Setup
if ~exist('eeg_checkset', 'file')
    addpath('~/eeglab');
end
eeglab nogui;

% Add BrainBeats functions to path (clean_rr lives there)
bbPath = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(bbPath, 'functions'));

testDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(testDir, 'data');
figDir  = fullfile(testDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% Load beat detections
beats = readtable(fullfile(dataDir, 'test_beat_detections.csv'));
beat_idx = beats.hrBeat == 1;
beat_sec = beats.timestamp_sec(beat_idx);
fprintf('Pre-detected beats: %d (mean IBI=%.0f ms)\n', length(beat_sec), mean(diff(beat_sec))*1000);

%% Convert to RR intervals (same as brainbeats_process rr branch)
target_srate = 250;

% Load EEG to get recording length for sample conversion
eeg_raw = readtable(fullfile(dataDir, 'test_eeg.csv'));
[eeg_ts, ~] = sort(eeg_raw.ts);
[eeg_ts, ~] = unique(eeg_ts);
t_start = eeg_ts(1);
t_end   = eeg_ts(end);
npnts = length(t_start:(1/target_srate):t_end);

Rpeaks = round(beat_sec * target_srate);
Rpeaks(Rpeaks < 1) = [];
Rpeaks(Rpeaks > npnts) = [];

RR   = diff(Rpeaks) / target_srate;       % seconds
RR_t = Rpeaks(1:end-1) / target_srate;    % seconds
fprintf('RR intervals: mean=%.0f ms, std=%.0f ms, n=%d\n', mean(RR)*1000, std(RR)*1000, length(RR));

%% Run clean_rr
sig_dummy = zeros(1, npnts);
[NN, NN_t, idx_rem, idx_interp] = clean_rr(RR_t, RR, sig_dummy(Rpeaks(1:end-1))');
fprintf('clean_rr: removed %d, interpolated %d out of %d intervals\n', sum(idx_rem), sum(idx_interp), length(RR));

% Reconstruct post-clean_rr peaks (same logic as brainbeats_process)
Rpeaks = Rpeaks(:)';  % ensure row vector
Rpeaks_clean = Rpeaks(2:end);
Rpeaks_clean(idx_rem) = [];
Rpeaks_after = [Rpeaks(1) Rpeaks_clean];

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
