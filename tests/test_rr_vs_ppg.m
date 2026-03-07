% test_rr_vs_ppg.m
% Compare heartbeat detection: BrainBeats PPG peak detection vs pre-detected beats.
%
% Loads Muse 2 EEG + optics data, runs BrainBeats in two modes:
%   1. 'ppg' mode — detect peaks from raw optical PPG signal
%   2. 'rr'  mode — use pre-detected beat latencies from beat_detections.csv
% Then overlays the detected peaks and HEP waveforms for comparison.

clear; close all;

%% Setup
if ~exist('eeg_checkset', 'file')
    addpath('~/eeglab');
end
eeglab nogui;

testDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(testDir, 'data');
figDir  = fullfile(testDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% Load raw data
fprintf('=== Loading raw data ===\n');
eeg_raw = readtable(fullfile(dataDir, 'test_eeg.csv'));
opt_raw = readtable(fullfile(dataDir, 'test_optics.csv'));
beats   = readtable(fullfile(dataDir, 'test_beat_detections.csv'));

%% Build EEGLAB dataset with EEG + PPG
target_srate = 250;

% EEG timestamps: Muse sends in packets, timestamps may be non-monotonic.
% Sort and remove duplicates to ensure interp1 works.
[eeg_ts, eeg_order] = sort(eeg_raw.ts);
[eeg_ts, uidx] = unique(eeg_ts);
eeg_order = eeg_order(uidx);

% Same for optics
[opt_ts, opt_order] = sort(opt_raw.ts);
[opt_ts, uidx_o] = unique(opt_ts);
opt_order = opt_order(uidx_o);

% Uniform time grid
t_start = eeg_ts(1);
t_end   = eeg_ts(end);
uniform_time = t_start:(1/target_srate):t_end;

% Resample EEG (4 channels) to uniform grid
eeg_ch = zeros(4, length(uniform_time));
for ch = 1:4
    col = sprintf('ch%d', ch);
    eeg_ch(ch,:) = interp1(eeg_ts, eeg_raw.(col)(eeg_order), uniform_time, 'pchip', 'extrap');
end

% Resample optics ch1 (PPG) to same uniform grid
ppg_ch = interp1(opt_ts, opt_raw.ch1(opt_order), uniform_time, 'pchip', 'extrap');

% Build EEG-only dataset (4 channels, no PPG) for approach A
EEG_only = pop_importdata('data', eeg_ch, 'srate', target_srate, 'setname', 'test_eeg_only');
eeg_labels = {'TP9', 'AF7', 'AF8', 'TP10'};
for i = 1:4
    EEG_only.chanlocs(i).labels = eeg_labels{i};
end
EEG_only = eeg_checkset(EEG_only);

% Build EEG+PPG dataset (5 channels) for approach B
all_data = [eeg_ch; ppg_ch];
EEG_with_ppg = pop_importdata('data', all_data, 'srate', target_srate, 'setname', 'test_eeg_ppg');
labels = {'TP9', 'AF7', 'AF8', 'TP10', 'PPG'};
for i = 1:length(labels)
    EEG_with_ppg.chanlocs(i).labels = labels{i};
end
EEG_with_ppg = eeg_checkset(EEG_with_ppg);
fprintf('Dataset: %d ch (+PPG), %d samples, %.1f sec, %d Hz\n', ...
    EEG_with_ppg.nbchan, EEG_with_ppg.pnts, EEG_with_ppg.pnts/EEG_with_ppg.srate, EEG_with_ppg.srate);

%% Extract beat latencies from beat_detections.csv
% timestamp_sec is relative time from recording start (same timebase as optics)
beat_idx = beats.hrBeat == 1;
beat_sec = beats.timestamp_sec(beat_idx);
fprintf('Pre-detected beats: %d (mean IBI=%.0f ms)\n', length(beat_sec), mean(diff(beat_sec))*1000);

%% Pre-filter EEG (BrainBeats clean_eeg requires channel locations for ICA/ICLabel;
%  Muse 4-ch data lacks them, so we filter manually and disable clean_eeg)
fprintf('Pre-filtering EEG 1-40 Hz...\n');
EEG_only = pop_eegfiltnew(EEG_only, 'locutoff', 1, 'hicutoff', 40);
EEG_with_ppg = pop_eegfiltnew(EEG_with_ppg, 'locutoff', 1, 'hicutoff', 40, 'channels', 1:4);

%% Approach A: Run BrainBeats with pre-detected beats ('rr' mode)
fprintf('\n=== Approach A: pre-detected beats (rr mode) ===\n');
EEG_A = brainbeats_process(EEG_only, ...
    'analysis', 'hep', ...
    'heart_signal', 'rr', ...
    'beat_latencies', beat_sec, ...
    'clean_eeg', false, ...
    'vis_cleaning', false, ...
    'vis_outputs', false, ...
    'save', false);

% Extract R-peak latencies from events
rpeaks_A_idx = strcmp({EEG_A.event.type}, 'R-peak');
rpeaks_A_lat = [EEG_A.event(rpeaks_A_idx).latency] / target_srate;  % seconds
fprintf('Approach A: %d R-peak events, %d epochs\n', sum(rpeaks_A_idx), EEG_A.trials);

%% Approach B: Run BrainBeats with raw PPG ('ppg' mode)
fprintf('\n=== Approach B: raw optical PPG (ppg mode) ===\n');
ppg_ok = true;
try
    EEG_B = brainbeats_process(EEG_with_ppg, ...
        'analysis', 'hep', ...
        'heart_signal', 'ppg', ...
        'heart_channels', {'PPG'}, ...
        'clean_eeg', false, ...
        'vis_cleaning', false, ...
        'vis_outputs', false, ...
        'save', false);

    % Extract R-peak latencies from events
    rpeaks_B_idx = strcmp({EEG_B.event.type}, 'R-peak');
    rpeaks_B_lat = [EEG_B.event(rpeaks_B_idx).latency] / target_srate;
    fprintf('Approach B: %d R-peak events, %d epochs\n', sum(rpeaks_B_idx), EEG_B.trials);
catch ME
    ppg_ok = false;
    fprintf('Approach B FAILED: %s\n', ME.message);
    fprintf('(Muse optics signal may not be compatible with BrainBeats PPG peak detector)\n');
    rpeaks_B_lat = [];
end

%% Figure 1: Peaks on PPG signal (first 30 seconds)
fig1 = figure('Position', [100 100 1200 500], 'Color', 'w');
t_ppg = (0:length(ppg_ch)-1) / target_srate;
xlims = [5 35];
idx_show = t_ppg >= xlims(1) & t_ppg <= xlims(2);

plot(t_ppg(idx_show), ppg_ch(idx_show), 'k', 'LineWidth', 0.5); hold on;

% A-peaks (from BrainBeats rr mode, after clean_rr)
in_A = rpeaks_A_lat >= xlims(1) & rpeaks_A_lat <= xlims(2);
rpA_samp = round(rpeaks_A_lat(in_A) * target_srate);
rpA_samp = max(1, min(length(ppg_ch), rpA_samp));
stem(rpeaks_A_lat(in_A), ppg_ch(rpA_samp), 'r', 'filled', 'MarkerSize', 6, 'LineWidth', 1);

% Original raw beat_detections (before clean_rr)
in_raw = beat_sec >= xlims(1) & beat_sec <= xlims(2);
raw_samp = round(beat_sec(in_raw) * target_srate);
raw_samp = max(1, min(length(ppg_ch), raw_samp));
plot(beat_sec(in_raw), ppg_ch(raw_samp), 'gv', 'MarkerSize', 5, 'MarkerFaceColor', 'g');

lgd = {'PPG signal', 'A: rr mode (after clean\_rr)', 'Raw beat\_detections'};
if ppg_ok
    in_B = rpeaks_B_lat >= xlims(1) & rpeaks_B_lat <= xlims(2);
    rpB_samp = round(rpeaks_B_lat(in_B) * target_srate);
    rpB_samp = max(1, min(length(ppg_ch), rpB_samp));
    stem(rpeaks_B_lat(in_B), ppg_ch(rpB_samp) * 0.98, 'b', 'filled', 'MarkerSize', 5);
    lgd{end+1} = 'B: PPG-detected';
end
xlabel('Time (s)'); ylabel('PPG (a.u.)');
title('Peak Detection: pre-detected beats vs PPG signal');
legend(lgd, 'Location', 'best');
xlim(xlims);
set(gca, 'FontSize', 11);
saveas(fig1, fullfile(figDir, 'fig1_peak_overlay.png'));
fprintf('Saved: fig1_peak_overlay.png\n');

%% Figure 2: IBI time series
fig2 = figure('Position', [100 100 1200 400], 'Color', 'w');
ibi_A = diff(rpeaks_A_lat) * 1000;
ibi_raw = diff(beat_sec) * 1000;

plot(beat_sec(2:end), ibi_raw, 'g.-', 'MarkerSize', 3); hold on;
plot(rpeaks_A_lat(2:end), ibi_A, 'r.-', 'MarkerSize', 4);
lgd2 = {'Raw beat\_detections', 'A: rr mode (after clean\_rr)'};
if ppg_ok
    ibi_B = diff(rpeaks_B_lat) * 1000;
    plot(rpeaks_B_lat(2:end), ibi_B, 'b.-', 'MarkerSize', 4);
    lgd2{end+1} = 'B: PPG-detected';
end
xlabel('Time (s)'); ylabel('IBI (ms)');
title('Inter-Beat Interval Time Series');
legend(lgd2, 'Location', 'best');
ylim([400 1200]);
set(gca, 'FontSize', 11);
saveas(fig2, fullfile(figDir, 'fig2_ibi_timeseries.png'));
fprintf('Saved: fig2_ibi_timeseries.png\n');

%% Figure 3: HEP waveforms from approach A
fig3 = figure('Position', [100 100 1200 800], 'Color', 'w');
chanNames = {'TP9', 'AF7', 'AF8', 'TP10'};
nChan = min(4, size(EEG_A.data,1));

for ch = 1:nChan
    subplot(2, 2, ch);
    hep_A = mean(EEG_A.data(ch,:,:), 3);
    plot(EEG_A.times, hep_A, 'r', 'LineWidth', 1.5); hold on;
    if ppg_ok
        hep_B = mean(EEG_B.data(ch,:,:), 3);
        plot(EEG_B.times, hep_B, 'b', 'LineWidth', 1.5);
    end
    xline(0, 'k--');
    yl = ylim;
    patch([200 350 350 200], [yl(1) yl(1) yl(2) yl(2)], ...
        [1 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([400 600 600 400], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.9 0.9 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    xlabel('Time (ms)'); ylabel('\muV');
    title(chanNames{ch});
    if ch == 1
        if ppg_ok
            legend({'A: pre-detected', 'B: PPG', 'R-peak', 'Early', 'Late'}, ...
                'Location', 'best', 'FontSize', 8);
        else
            legend({'A: pre-detected', 'R-peak', 'Early (200-350)', 'Late (400-600)'}, ...
                'Location', 'best', 'FontSize', 8);
        end
    end
    set(gca, 'FontSize', 10);
end
if ppg_ok
    sgtitle(sprintf('HEP: A (%d epochs) vs B (%d epochs)', EEG_A.trials, EEG_B.trials), 'FontSize', 13);
else
    sgtitle(sprintf('HEP from pre-detected beats (%d epochs)', EEG_A.trials), 'FontSize', 13);
end
saveas(fig3, fullfile(figDir, 'fig3_hep_waveforms.png'));
fprintf('Saved: fig3_hep_waveforms.png\n');

%% Figure 4: IBI distribution
fig4 = figure('Position', [100 100 600 400], 'Color', 'w');
histogram(ibi_raw, 40, 'FaceColor', [0.6 0.9 0.6], 'FaceAlpha', 0.5); hold on;
histogram(ibi_A, 40, 'FaceColor', [0.9 0.4 0.4], 'FaceAlpha', 0.5);
xlabel('IBI (ms)'); ylabel('Count');
title('IBI Distribution');
legend({'Raw beat\_detections', 'After clean\_rr'}, 'Location', 'best');
set(gca, 'FontSize', 11);
saveas(fig4, fullfile(figDir, 'fig4_ibi_distribution.png'));
fprintf('Saved: fig4_ibi_distribution.png\n');

%% Summary
n_A = length(rpeaks_A_lat);
fprintf('\n=== SUMMARY ===\n');
fprintf('Raw beat_detections: %d beats\n', length(beat_sec));
fprintf('Approach A (rr mode, after clean_rr): %d peaks, %d HEP epochs\n', n_A, EEG_A.trials);
if ppg_ok
    n_B = length(rpeaks_B_lat);
    fprintf('Approach B (ppg mode): %d peaks, %d HEP epochs\n', n_B, EEG_B.trials);
else
    fprintf('Approach B (ppg mode): FAILED (Muse optics incompatible with BrainBeats PPG detector)\n');
end

% HEP amplitude summary
fprintf('\nHEP amplitudes (Approach A):\n');
for ch = 1:nChan
    early_A = mean(EEG_A.data(ch, EEG_A.times>=200 & EEG_A.times<=350, :), [2 3]);
    late_A  = mean(EEG_A.data(ch, EEG_A.times>=400 & EEG_A.times<=600, :), [2 3]);
    fprintf('  %s: early=%.2f uV, late=%.2f uV\n', chanNames{ch}, early_A, late_A);
end

fprintf('\nAll figures saved to: %s\n', figDir);
fprintf('Test complete.\n');
