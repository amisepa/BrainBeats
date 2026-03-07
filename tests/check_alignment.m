% check_alignment.m — Investigate timestamp alignment between EEG, optics, and beat_detections

eeg = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'test_eeg.csv'));
opt = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'test_optics.csv'));
beats = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'test_beat_detections.csv'));

fprintf('=== Row counts ===\n');
fprintf('EEG rows: %d, Optics rows: %d, Beats rows: %d\n', height(eeg), height(opt), height(beats));

fprintf('\n=== Unix timestamp ranges ===\n');
fprintf('EEG ts:     [%.6f .. %.6f]\n', eeg.ts(1), eeg.ts(end));
fprintf('Optics ts:  [%.6f .. %.6f]\n', opt.ts(1), opt.ts(end));
fprintf('Optics start - EEG start = %.4f sec\n', opt.ts(1) - eeg.ts(1));
fprintf('Optics end   - EEG end   = %.4f sec\n', opt.ts(end) - eeg.ts(end));

fprintf('\n=== EEG timestamp intervals (diff) ===\n');
eeg_dt = diff(eeg.ts);
fprintf('mean=%.6f s, std=%.6f s, min=%.6f s, max=%.6f s\n', mean(eeg_dt), std(eeg_dt), min(eeg_dt), max(eeg_dt));
fprintf('Negative diffs: %d, Zero diffs: %d\n', sum(eeg_dt < 0), sum(eeg_dt == 0));
fprintf('Expected dt at 256Hz = %.6f s\n', 1/256);

fprintf('\n=== Optics timestamp intervals (diff) ===\n');
opt_dt = diff(opt.ts);
fprintf('mean=%.6f s, std=%.6f s, min=%.6f s, max=%.6f s\n', mean(opt_dt), std(opt_dt), min(opt_dt), max(opt_dt));
fprintf('Negative diffs: %d, Zero diffs: %d\n', sum(opt_dt < 0), sum(opt_dt == 0));
fprintf('Expected dt at 64Hz = %.6f s\n', 1/64);

fprintf('\n=== Beat detections timestamp_sec ===\n');
fprintf('Range: [%.6f .. %.6f]\n', beats.timestamp_sec(1), beats.timestamp_sec(end));
beat_dt = diff(beats.timestamp_sec);
fprintf('Interval mean=%.6f s, std=%.6f s, min=%.6f s, max=%.6f s\n', mean(beat_dt), std(beat_dt), min(beat_dt), max(beat_dt));
fprintf('Negative diffs: %d\n', sum(beat_dt < 0));

fprintf('\n=== Row-level alignment: beats.timestamp_sec vs optics relative time ===\n');
opt_relative = opt.ts - opt.ts(1);
fprintf('beats.timestamp_sec(1)   = %.6f,  opt_relative(1)   = %.6f\n', beats.timestamp_sec(1), opt_relative(1));
fprintf('beats.timestamp_sec(end) = %.6f,  opt_relative(end) = %.6f\n', beats.timestamp_sec(end), opt_relative(end));
if height(beats) == height(opt)
    max_diff = max(abs(beats.timestamp_sec - opt_relative));
    fprintf('Max |beats.timestamp_sec - opt_relative| = %.6f sec\n', max_diff);
else
    fprintf('Row counts differ — cannot do row-level comparison\n');
end

fprintf('\n=== Critical offset for test script ===\n');
offset_sec = opt.ts(1) - eeg.ts(1);
fprintf('beat_detections t=0 is optics start.\n');
fprintf('Optics start is %.4f sec AFTER EEG start.\n', offset_sec);
fprintf('In the test, t_ppg=0 maps to EEG start (uniform_time begins at eeg_ts(1)).\n');
fprintf('So beat_sec should be shifted by +%.4f sec to align with the PPG plot.\n', offset_sec);

% Show first 10 beat times and where they fall on PPG
fprintf('\n=== First 10 beat times (raw vs corrected) ===\n');
bi = beats.hrBeat == 1;
bsec = beats.timestamp_sec(bi);
fprintf('%8s  %12s  %12s\n', 'Beat#', 'Raw(s)', 'Corrected(s)');
for k = 1:min(10, length(bsec))
    fprintf('%8d  %12.4f  %12.4f\n', k, bsec(k), bsec(k) + offset_sec);
end

fprintf('\nDone.\n');
