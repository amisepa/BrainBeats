function compare_badepoch_metrics(files)
% COMPARE_BADEPOCH_METRICS - Bad-epoch detection for HEP (FIND_BADTRIALS):
% the current metrics (amplitude, and high-frequency noise = data minus their
% 45-Hz low-pass, on the padded epochs of the 0.5-30 Hz data) against the
% same on the analysis window only, and the high-frequency metric on the
% broadband copy the ICA is fitted on (1-Hz high-pass, line noise notched).
%
% Muscle-like bursts (300 ms of 30-100 Hz noise, 15 uV RMS, on 10 temporal
% channels) are added to the raw EEG at 15 random heartbeats (centered 0 to
% 300 ms after the R-peak). Hits: flagged epochs whose analysis window (-300
% to 600 ms) holds a burst; collateral: flagged epochs whose window holds
% none but whose padded epoch (650 ms on each side) does; also the epochs
% flagged without any burst. isoutlier 'grubbs' as in CLEAN_EEG. R-peaks
% from the ECG (BrainBeats HRV run), average reference, no IBI rejection.
%
% Usage:
%   compare_badepoch_metrics({'sample_data/sub-047_rest.set'})
%
% Results (2026-11-05; sub-047, dataset.set, ds003838 sub-062, 033, 040,
% 070): bursts caught (of 15) / other epochs flagged (no burst nearby):
%   v1.6 metrics (padded epochs, MAD of the >45 Hz residual of the 0.5-30
%   Hz data): 0, 0, 2, 0, 0, 0 / 4, 16, 14, 2, 0, 2. The MAD over time
%   ignores a 300-ms burst, and the 0.5-30 Hz data hold nothing above 45 Hz.
%   FIND_BADTRIALS now (analysis window; 90th percentile over channels of
%   the RMS of the >45 Hz residual of the broadband copy): 15, 15, 15, 1,
%   15, 15 / 3, 17, 12, 8, 1, 8 (sub-033 and sub-070 have much muscle
%   activity). Maximum over channels or a log transform did worse (0/15 on
%   3 recordings).
%
% Copyright (C) - Cedric Cannard, 2026

if ~exist('pop_loadset','file'), eeglab; close; end
bb = fileparts(fileparts(mfilename('fullpath')));
addpath(bb, fullfile(bb,'functions'));
temporal = {'T7' 'T8' 'TP7' 'TP8' 'FT7' 'FT8' 'F7' 'F8' 'P7' 'P8'};
files = cellstr(files);
for k = 1:numel(files)
    [fp, fn, ext] = fileparts(files{k});
    EEG = pop_loadset('filename',[fn ext],'filepath',fp);
    fs = EEG.srate;
    E = brainbeats_process(pop_select(EEG,'channel',{'ECG'}),'analysis','features','eeg','off', ...
        'heart_signal','ECG','heart_channels',{'ECG'},'hrv_features',{'time'},'save',0, ...
        'vis_cleaning',0,'vis_outputs',0,'gong',0);
    NN = E.brainbeats.preprocessings.NN(:);
    beats = round(([0; cumsum(NN)] + E.brainbeats.preprocessings.NN_times(1) - NN(1)) * fs) + 1;
    pad = round(0.65*fs); w0 = round(-0.3*fs); w1 = round(0.6*fs) - 1;
    beats = beats(beats + w0 - pad > 1 & beats + w1 + pad < EEG.pnts);
    X = pop_select(EEG,'nochannel',intersect({EEG.chanlocs.labels}, {'ECG' 'PPG'}));
    for withBursts = [false true]
        Y = X;
        rng(k); tb = [];
        if withBursts
            ch = ismember({Y.chanlocs.labels}, temporal);
            pick = sort(randperm(numel(beats), 15));
            [bh, ah] = butter(4, [30 min(100, 0.45*fs)]/(fs/2));
            n = round(0.3*fs);
            for b = beats(pick)'
                c = b + round(0.3*rand*fs);                 % burst center 0-300 ms after the R-peak
                s = filtfilt(bh, ah, randn(sum(ch), n + 200)')';
                s = s(:, 101:100+n); s = 15 * s ./ std(s, 0, 2) .* hann(n)' / sqrt(mean(hann(n).^2));
                Y.data(ch, c - floor(n/2) + (0:n-1)) = Y.data(ch, c - floor(n/2) + (0:n-1)) + s;
                tb(end+1) = c; %#ok<AGROW>
            end
        end
        D = pop_reref(pop_eegfiltnew(pop_eegfiltnew(Y,'locutoff',0.5),'hicutoff',30), []);      % HEP data
        C = pop_reref(pop_eegfiltnew(Y,'locutoff',1), []);                                    % ICA copy
        for f = 50:50:fs/2-3.5
            C = pop_eegfiltnew(C,'locutoff',f-3,'hicutoff',f+3,'revfilt',1,'filtorder',2*round(fs));
        end
        idx = (w0 - pad):(w1 + pad); inWin = idx >= w0 & idx <= w1;
        ep = @(Z) reshape(double(Z.data(:, beats + idx)), Z.nbchan, numel(beats), []);   % ch x beats x times
        De = permute(ep(D), [1 3 2]); Ce = permute(ep(C), [1 3 2]);                       % ch x times x beats
        M = struct();
        [M.amp_padded, M.hf_padded] = metrics(De, fs);                 % current
        [M.amp_window, M.hf_window] = metrics(De(:,inWin,:), fs);
        [~, M.hf_copy_window, M.hfrms_copy_max, M.hfrms_copy_p90] = metrics(Ce(:,inWin,:), fs);
        M.log_amp_window = log(M.amp_window);
        M.log_hfrms_copy_max = log(M.hfrms_copy_max);
        M.log_hfrms_copy_p90 = log(M.hfrms_copy_p90);
        names = fieldnames(M);
        % FIND_BADTRIALS as implemented (high-frequency metric on the copy,
        % the whole padded epoch filtered, then the analysis window)
        Dep = struct('data', De, 'times', idx/fs*1000, 'srate', fs, 'pnts', numel(idx), 'trials', numel(beats));
        nowBad = false(1, numel(beats));
        nowBad(find_badTrials(Dep, 'grubbs', false, struct('win', [-300 600], 'hf', Ce, 'lowpass', 30))) = true;
        if ~withBursts
            base = structfun(@(m) sum(isoutlier(m,'grubbs')), M);
            baseNow = sum(nowBad);
            continue
        end
        tgt = false(1, numel(beats)); col = tgt;
        for t = tb
            tgt = tgt | (t >= beats' + w0 & t <= beats' + w1);
            col = col | (t >= beats' + w0 - pad & t <= beats' + w1 + pad);
        end
        col = col & ~tgt;
        fprintf('BADEP %s: %d heartbeats, 15 bursts in %d windows\n', fn, numel(beats), sum(tgt));
        for i = 1:numel(names)
            bad = isoutlier(M.(names{i}), 'grubbs');
            fprintf('BADEP   %-15s flagged %3d (no bursts: %3d): hits %2d/%d, collateral %2d, other %2d\n', names{i}, ...
                sum(bad), base(i), sum(bad & tgt), sum(tgt), sum(bad & col), sum(bad & ~tgt & ~col));
        end
        fprintf('BADEP   find_badTrials now flags %d epochs without bursts\n', baseNow);
        rules = {'v1.6 (amp_padded | hf_padded)', isoutlier(M.amp_padded,'grubbs') | isoutlier(M.hf_padded,'grubbs')
                 'amp_window | hfrms_copy_p90 (window filtered)', isoutlier(M.amp_window,'grubbs') | isoutlier(M.hfrms_copy_p90,'grubbs')
                 'find_badTrials now', nowBad};
        for i = 1:size(rules,1)
            bad = rules{i,2};
            fprintf('BADEP   rule %-46s hits %2d/%d, collateral %2d, other %2d\n', rules{i,1}, sum(bad & tgt), sum(tgt), ...
                sum(bad & col), sum(bad & ~tgt & ~col));
        end
    end
end
end

function [amp, hf, hfmax, hfp90] = metrics(X, fs)
% FIND_BADTRIALS metrics of epochs X (channels x times x epochs): amplitude,
% high-frequency noise (MAD over time, RMS over channels; current), and the
% RMS over time of the high-frequency residual, maximum and 90th percentile
% over channels
b = design_fir(100, [2*[0 45 50]/fs 1], [1 1 0 0]);
n = size(X,3); amp = nan(1,n); hf = amp; hfmax = amp; hfp90 = amp;
for e = 1:n
    x = X(:,:,e);
    amp(e) = rms(rms(x, 2));
    r = x - filtfilt_fast(b, 1, x')';
    hf(e) = rms(mad(r, 0, 2));
    q = rms(r, 2);
    hfmax(e) = max(q); hfp90(e) = prctile(q, 90);
end
end
