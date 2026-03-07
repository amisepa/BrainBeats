% test_rr_vs_ppg.m
% Compare heartbeat detection: BrainBeats PPG peak detection vs pre-detected beats.
%
% Loads Muse 2 EEG + optics data, runs BrainBeats in two modes:
%   1. 'rr'  mode — use pre-detected beat latencies from beat_detections.csv
%   2. 'ppg' mode — detect peaks from raw optical PPG signal
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

% Use corrected timestamps from beat_detections as the optics time axis.
% beat_detections has the same row count as optics (28152) and uses a
% corrected ideal 64 Hz clock, while raw optics ts has Bluetooth jitter
% (28% out-of-order, up to 1.8 s cumulative drift).
opt_time_corrected = beats.timestamp_sec;  % ideal 64 Hz grid, relative seconds
fprintf('Optics time axis: using corrected timestamps from beat_detections (ideal 64 Hz)\n');

% EEG timestamps: Muse sends in packets, timestamps may be non-monotonic.
% Sort and remove duplicates, then convert to relative seconds.
[eeg_ts, eeg_order] = sort(eeg_raw.ts);
[eeg_ts, uidx] = unique(eeg_ts);
eeg_order = eeg_order(uidx);
eeg_ts_rel = eeg_ts - eeg_ts(1);  % relative seconds from EEG start

% Optics starts ~36 ms after EEG; shift optics time axis to EEG-relative.
% This aligns both streams to a common t=0 (EEG start).
eeg_opt_offset = opt_raw.ts(1) - eeg_raw.ts(1);  % ~0.036 s
opt_time_aligned = opt_time_corrected + eeg_opt_offset;

% Uniform time grid in relative seconds
t_start = 0;
t_end   = eeg_ts_rel(end);
uniform_time = t_start:(1/target_srate):t_end;

% Resample EEG (4 channels) to uniform grid
eeg_ch = zeros(4, length(uniform_time));
for ch = 1:4
    col = sprintf('ch%d', ch);
    eeg_ch(ch,:) = interp1(eeg_ts_rel, eeg_raw.(col)(eeg_order), uniform_time, 'pchip', 'extrap');
end

% Resample optics ch1 (PPG) to same uniform grid using corrected timestamps
ppg_ch = interp1(opt_time_aligned, opt_raw.ch1, uniform_time, 'pchip', 'extrap');

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
% Use corrected timestamps, shifted to EEG-relative time (same as PPG)
beat_idx = beats.hrBeat == 1;
beat_sec = beats.timestamp_sec(beat_idx) + eeg_opt_offset;
fprintf('Pre-detected beats: %d (mean IBI=%.0f ms)\n', length(beat_sec), mean(diff(beat_sec))*1000);

%% Pre-filter EEG (BrainBeats clean_eeg requires channel locations for ICA/ICLabel;
%  Muse 4-ch data lacks them, so we filter manually and disable clean_eeg)
fprintf('Pre-filtering EEG 1-40 Hz...\n');
EEG_only = pop_eegfiltnew(EEG_only, 'locutoff', 1, 'hicutoff', 40);
EEG_with_ppg = pop_eegfiltnew(EEG_with_ppg, 'locutoff', 1, 'hicutoff', 40, 'channels', 1:4);

%% Approach A: Run BrainBeats with pre-detected beats ('rr' mode, no clean_rr)
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

%% Approach C: smoothed findpeaks on PPG signal
fprintf('\n=== Approach C: findpeaks on smoothed PPG ===\n');
t_ppg = (0:length(ppg_ch)-1) / target_srate;

% Ground truth: pre-detected beats in 21-35s window
gt_mask = beat_sec >= 21 & beat_sec <= 35;
gt_sec  = beat_sec(gt_mask);
fprintf('Ground truth window [21-35s]: %d beats\n', length(gt_sec));

% Grid search over smoothing window and MinPeakProminence
smooth_wins  = [3 5 7 9 11 15 19 25 31];       % movmean window (samples)
prominences  = [0.0005 0.001 0.002 0.003 0.005 0.007 0.01 0.015 0.02];
min_dist_samp = round(0.3 * target_srate);      % 75 samples = 0.3 s
tol_sec = 0.35;                                  % 350 ms tolerance (PPG peak vs beat onset phase diff)

best_f1 = 0; best_sw = 3; best_prom = 0.001;    % safe defaults
for sw = smooth_wins
    ppg_smooth = movmean(ppg_ch, sw);
    for prom = prominences
        [~, locs] = findpeaks(-ppg_smooth, 'MinPeakDistance', min_dist_samp, ...
            'MinPeakProminence', prom);
        det_sec = (locs - 1) / target_srate;

        % Match detected peaks to ground truth within tolerance
        det_gt = det_sec(det_sec >= 21 & det_sec <= 35);
        hits = 0;
        gt_matched = false(size(gt_sec));
        for k = 1:length(det_gt)
            diffs = abs(gt_sec - det_gt(k));
            [mn, mi] = min(diffs);
            if mn <= tol_sec && ~gt_matched(mi)
                hits = hits + 1;
                gt_matched(mi) = true;
            end
        end
        prec = hits / max(length(det_gt), 1);
        rec  = hits / max(length(gt_sec), 1);
        f1   = 2 * prec * rec / max(prec + rec, eps);
        if f1 > best_f1
            best_f1 = f1; best_sw = sw; best_prom = prom;
        end
    end
end
fprintf('Best params: smoothing=%d samples, MinPeakProminence=%.4f (F1=%.3f)\n', ...
    best_sw, best_prom, best_f1);

% Apply best parameters to full signal
ppg_smooth_best = movmean(ppg_ch, best_sw);
[fp_vals, fp_locs] = findpeaks(-ppg_smooth_best, 'MinPeakDistance', min_dist_samp*1.25, ...
    'MinPeakProminence', best_prom);
fp_sec = (fp_locs - 1) / target_srate;
fprintf('Approach C: %d peaks detected over full recording\n', length(fp_sec));

%% Figure 1: Peaks on PPG signal (first 30 seconds) — three panels
fig1 = figure('Position', [100 100 1200 900], 'Color', 'w');
xlims = [5 35];
idx_show = t_ppg >= xlims(1) & t_ppg <= xlims(2);

% Panel 1: pre-detected beats (green)
subplot(3,1,1);
plot(t_ppg(idx_show), ppg_ch(idx_show), 'k', 'LineWidth', 0.5); hold on;
in_raw = beat_sec >= xlims(1) & beat_sec <= xlims(2);
raw_samp = round(beat_sec(in_raw) * target_srate);
raw_samp = max(1, min(length(ppg_ch), raw_samp));
stem(beat_sec(in_raw), ppg_ch(raw_samp), 'g', 'filled', 'MarkerSize', 5, 'LineWidth', 1);
ylabel('PPG (a.u.)');
title('Pre-detected beats (from beat\_detections.csv)');
legend({'PPG signal', 'Pre-detected beats'}, 'Location', 'best');
xlim(xlims); ylim([1.06 1.12]);
set(gca, 'FontSize', 11);

% Panel 2: BrainBeats rr-mode peaks (red)
subplot(3,1,2);
plot(t_ppg(idx_show), ppg_ch(idx_show), 'k', 'LineWidth', 0.5); hold on;
in_A = rpeaks_A_lat >= xlims(1) & rpeaks_A_lat <= xlims(2);
rpA_samp = round(rpeaks_A_lat(in_A) * target_srate);
rpA_samp = max(1, min(length(ppg_ch), rpA_samp));
stem(rpeaks_A_lat(in_A), ppg_ch(rpA_samp), 'r', 'filled', 'MarkerSize', 5, 'LineWidth', 1);
ylabel('PPG (a.u.)');
title('BrainBeats rr mode (after run\_HEP)');
legend({'PPG signal', 'BrainBeats R-peak events'}, 'Location', 'best');
xlim(xlims); ylim([1.06 1.12]);
set(gca, 'FontSize', 11);

% Panel 3: findpeaks on smoothed PPG (blue)
subplot(3,1,3);
plot(t_ppg(idx_show), ppg_ch(idx_show), 'k', 'LineWidth', 0.5); hold on;
plot(t_ppg(idx_show), ppg_smooth_best(idx_show), 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
in_fp = fp_sec >= xlims(1) & fp_sec <= xlims(2);
stem(fp_sec(in_fp), ppg_ch(fp_locs(in_fp)), 'b', 'filled', 'MarkerSize', 5, 'LineWidth', 1);
xlabel('Time (s)'); ylabel('PPG (a.u.)');
title(sprintf('findpeaks (smooth=%d, prom=%.4f, minDist=0.3s, F1=%.2f)', ...
    best_sw, best_prom, best_f1));
legend({'PPG raw', 'PPG smoothed', 'findpeaks'}, 'Location', 'best');
xlim(xlims); ylim([1.06 1.12]);
set(gca, 'FontSize', 11);

saveas(fig1, fullfile(figDir, 'fig1_peak_overlay.png'));
fprintf('Saved: fig1_peak_overlay.png\n');

%% Figure 2: IBI time series
fig2 = figure('Position', [100 100 1200 400], 'Color', 'w');
ibi_A = diff(rpeaks_A_lat) * 1000;

plot(rpeaks_A_lat(2:end), ibi_A, 'r.-', 'MarkerSize', 4); hold on;
lgd2 = {'A: rr mode (pre-detected)'};
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

%% Summary
n_A = length(rpeaks_A_lat);
fprintf('\n=== SUMMARY ===\n');
fprintf('Raw beat_detections: %d beats\n', length(beat_sec));
fprintf('Approach A (rr mode): %d peaks, %d HEP epochs\n', n_A, EEG_A.trials);
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
