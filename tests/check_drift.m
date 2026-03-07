% check_drift.m — Visualize timestamp drift between ideal beat_detections grid and actual optics timestamps

close all;
opt = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'test_optics.csv'));
beats = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'test_beat_detections.csv'));

% Optics: sort and get monotonic time
[opt_ts_sorted, ord] = sort(opt.ts);
[opt_ts_unique, uidx] = unique(opt_ts_sorted);

% Ideal 64 Hz grid (what beat_detections assumes)
ideal_time = (0:height(opt)-1)' / 64;

% Actual optics relative time (from Unix timestamps)
opt_relative = opt.ts - opt.ts(1);

% Drift = ideal - actual
drift = beats.timestamp_sec - opt_relative;

figDir = fullfile(fileparts(mfilename('fullpath')), 'figures');

fig = figure('Position', [100 100 1200 800], 'Color', 'w');

subplot(3,1,1);
plot(opt_relative, drift, 'b.', 'MarkerSize', 1);
xlabel('Optics time (s)'); ylabel('Drift (s)');
title('beat\_detections ideal time - optics actual time');
set(gca, 'FontSize', 11);

subplot(3,1,2);
opt_dt = diff(opt.ts);
plot(opt_relative(2:end), opt_dt * 1000, '.', 'MarkerSize', 1);
ylabel('dt (ms)'); xlabel('Optics time (s)');
title('Optics inter-sample intervals (raw, unsorted)');
ylim([-20 200]);
yline(1000/64, 'r--', 'LineWidth', 1);
legend({'Actual dt', 'Expected 15.625 ms'}, 'Location', 'best');
set(gca, 'FontSize', 11);

subplot(3,1,3);
% Show sorted optics timestamps to see if monotonic after sort
opt_dt_sorted = diff(opt_ts_sorted);
plot((opt_ts_sorted(2:end) - opt_ts_sorted(1)), opt_dt_sorted * 1000, '.', 'MarkerSize', 1);
ylabel('dt (ms)'); xlabel('Time (s)');
title('Optics inter-sample intervals (after sorting)');
ylim([-2 100]);
yline(1000/64, 'r--', 'LineWidth', 1);
legend({'Sorted dt', 'Expected 15.625 ms'}, 'Location', 'best');
set(gca, 'FontSize', 11);

saveas(fig, fullfile(figDir, 'timestamp_drift.png'));
fprintf('Saved: timestamp_drift.png\n');

% Also check: can we remap beat times using the actual optics timestamps?
% Since beats has same row count as optics, beat at row i should correspond
% to optics sample i. So the real time of beat i = opt.ts(i) - opt.ts(1)
bi = beats.hrBeat == 1;
beat_rows = find(bi);
beat_sec_ideal = beats.timestamp_sec(bi);
beat_sec_actual = opt.ts(beat_rows) - opt.ts(1);

fprintf('\n=== Beat time comparison (ideal vs optics-based) ===\n');
fprintf('First beat: ideal=%.4f s, actual=%.4f s, diff=%.4f s\n', beat_sec_ideal(1), beat_sec_actual(1), beat_sec_ideal(1)-beat_sec_actual(1));
fprintf('Last beat:  ideal=%.4f s, actual=%.4f s, diff=%.4f s\n', beat_sec_ideal(end), beat_sec_actual(end), beat_sec_ideal(end)-beat_sec_actual(end));
fprintf('Max |diff|: %.4f s\n', max(abs(beat_sec_ideal - beat_sec_actual)));
fprintf('Mean |diff|: %.4f s\n', mean(abs(beat_sec_ideal - beat_sec_actual)));

fprintf('\nDone.\n');
