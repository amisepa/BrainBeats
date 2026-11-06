function check_surrogate_reference(files)
% CHECK_SURROGATE_REFERENCE - The HEP surrogate test with each reference of
% SIM_SURROGATE_REFERENCE (cycle average, pre-R -150 to -50 ms, pre-R -300
% to -200 ms, cycle without the QRS) on real recordings: the latencies where
% the ROI HEP differs from the surrogates (TFCE and max statistic, p < .05).
%
% Light preprocessing (EEG 0.5-30 Hz, average reference, no ICA, cardiac
% field left in), beats from the ECG (BrainBeats HRV run), heartbeats
% followed by the next one within 650 ms left out (as RUN_HEP with [-300 600]),
% frontocentral ROI, 1000 surrogate HEPs. Also prints the timing of the P and
% T waves of the ECG average and the share of heartbeats whose previous
% interval is below 750 ms (the previous T wave can then reach -300 ms).
%
% Usage:
%   check_surrogate_reference({'sample_data/sub-047_rest.set'})
%
% Results (2026-11-05; ds003838 recordings made with PREP_DS003838):
%   sub-047 (P wave -168 ms, previous IBI < 750 ms 0%): cycle, pre300 and
%     noQRS give the QRS, ~130-190, ~240-250 and ~400-430 ms; pre150 loses the
%     late positivity (TFCE) and adds 48-56 and 80-280 ms (its window holds
%     the P wave).
%   sub-062 (P wave -108 ms, previous IBI < 750 ms 82%): cycle 144-168,
%     248-280, 328-344, 380-464 ms and a pre-R point (-256 to -240); pre150
%     132-172, 244-280, 420-452 ms (no pre-R point); pre300 -52 to 16 and
%     44 to 308 ms (its window holds the previous heartbeat's late activity).
%   dataset.set (tutorial; P -132 ms, 57%): QRS only with all four.
%   sub-033, sub-070: none with all four; sub-040: QRS points at most.
% No fixed pre-R window is quiet in every recording: the P wave falls 100-170
% ms before the R-peak and the previous heartbeat reaches the epoch start at
% faster heart rates.
%
% Copyright (C) - Cedric Cannard, 2026

if ~exist('pop_loadset','file'), eeglab; close; end
bb = fileparts(fileparts(mfilename('fullpath')));
addpath(bb, fullfile(bb,'functions'));
files = cellstr(files);
roiLab = {'F1' 'F2' 'F3' 'F4' 'Fz' 'FC1' 'FC2' 'FC3' 'FC4' 'FCz' 'C1' 'C2' 'C3' 'C4' 'Cz'};
for k = 1:numel(files)
    [fp, fn, ext] = fileparts(files{k});
    EEG = pop_loadset('filename',[fn ext],'filepath',fp);
    fs = EEG.srate;
    E = brainbeats_process(pop_select(EEG,'channel',{'ECG'}),'analysis','features','eeg','off', ...
        'heart_signal','ECG','heart_channels',{'ECG'},'hrv_features',{'time'},'save',0, ...
        'vis_cleaning',0,'vis_outputs',0,'gong',0);
    NN = E.brainbeats.preprocessings.NN(:);
    beats = round(([0; cumsum(NN)] + E.brainbeats.preprocessings.NN_times(1) - NN(1)) * fs) + 1;
    ecg = pop_eegfiltnew(pop_select(EEG,'channel',{'ECG'}),'locutoff',0.5,'hicutoff',30);
    X = pop_select(EEG,'nochannel',intersect({EEG.chanlocs.labels}, {'ECG' 'PPG'}));
    X = pop_eegfiltnew(X,'locutoff',0.5,'hicutoff',30);
    X = pop_reref(X, []);
    roi = mean(double(X.data(ismember({X.chanlocs.labels}, roiLab),:)), 1);
    ibi = diff(beats) / fs * 1000;
    keep = [ibi >= 650; false];                       % next heartbeat after epoch end + 50 ms
    [tf, S] = compute_hep_tf(roi, fs, beats(keep), [-300 600], struct('tf',false,'nSurr',100, ...
        'nSurrHep',1000,'allBeats',beats,'keepNull',true));
    t = tf.hep_times;
    % ECG average: P (largest |deflection| at -300 to -60 ms) and T (150-500 ms)
    tb = tf.beats; e = zeros(size(t));
    for b = tb(:)', e = e + double(ecg.data(1, b + round(t/1000*fs))); end
    e = e / numel(tb); e = e - median(e(t >= -300 & t <= -60));
    [~, iP] = max(abs(e .* (t >= -300 & t <= -60))); [~, iT] = max(abs(e .* (t >= 150 & t <= 500)));
    kk = find(keep); kk = kk(kk > 1);
    prevShort = mean(ibi(kk-1) < 750);
    fprintf('\n%s: %d heartbeats (%d kept), P wave peak %g ms, T wave peak %g ms, previous IBI < 750 ms: %.0f%%\n', ...
        fn, numel(beats), numel(tb), t(iP), t(iT), 100*prevShort);
    refs = {'cycle', []; 'pre150', t >= -150 & t <= -50; 'pre300', t >= -300 & t <= -200; ...
        'noQRS', t < -50 | t > 100};
    for r = 1:size(refs,1)
        h = tf.hep; N = S.hep.null;
        if ~isempty(refs{r,2})
            h = h - mean(h(refs{r,2})); N = N - mean(N(:,refs{r,2},:), 2);
        end
        [pT, pM] = fwer(h, N);
        fprintf('  %-7s TFCE: %s | max statistic: %s\n', refs{r,1}, spans(t, pT < .05), spans(t, pM < .05));
    end
end
end

function s = spans(t, sig)
% Latency intervals (ms) where sig is true
d = diff([false sig false]); on = find(d == 1); off = find(d == -1) - 1;
if isempty(on), s = 'none'; return; end
s = strjoin(arrayfun(@(a,b) sprintf('%g to %g', t(a), t(b)), on, off, 'UniformOutput', false), ', ');
end

function [pT, pM] = fwer(h, N)
% TFCE and max-statistic p-values (as in SIM_SURROGATE_REFERENCE)
nS = size(N, 3);
mu = mean(N, 3); sd = std(N, 0, 3);
z = (h - mu) ./ sd;
d = reshape(N - mu, [], nS)';
SS = (nS-1) * sd.^2;
zn = (d * nS/(nS-1)) ./ sqrt(max(SS - d.^2 * nS/(nS-1), 0) / (nS-2));
maxZ = max(abs(zn), [], 2);
maxT = max(abs(tfce(zn)), [], 2);
pM = (sum(maxZ' >= abs(z(:)), 2) + 1)' / (nS + 1);
pT = (sum(maxT' >= abs(tfce(z))', 2) + 1)' / (nS + 1);
end

function T = tfce(Z)
T = zeros(size(Z));
for sgn = [1 -1]
    Zs = max(sgn*Z, 0);
    for h = 0.1:0.1:max(Zs(:))
        M = Zs >= h;
        c = cumsum(M, 2); L = c - cummax(c .* ~M, 2);
        c = cumsum(M, 2, 'reverse'); L = (L + c - cummax(c .* ~M, 2, 'reverse') - 1) .* M;
        T = T + sgn * h^2 * 0.1 * sqrt(L) .* M;
    end
end
end
