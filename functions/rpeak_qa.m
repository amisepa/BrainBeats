% RPEAK_QA - Check R-peak detections for T-wave or S-wave locking on one ECG
% recording. Report only: nothing is removed.
%
% A beat detected on the wrong deflection is still cardiac-timed, so it passes
% the SQI and the RR cleaning, but it jitters the HEP epoch alignment and
% smears the evoked response toward zero. Three checks on the detected beats:
%   1. Template correlation: each beat's ECG epoch against the median beat.
%      A subpopulation with low r = beats aligned to a different deflection.
%   2. R amplitude: fraction of beats whose peak-to-peak amplitude around
%      t = 0 is below half the median (small T-waves mixed with large R-waves).
%   3. RR alternation: lag-1 autocorrelation of the RR series. Strongly
%      negative = alternating short/long intervals = T-waves counted as beats
%      (normal RR has a weak positive lag-1 autocorrelation).
%
% Usage:
%   qa = rpeak_qa(ecg, peaks, fs)
%   qa = rpeak_qa(ecg, peaks, fs, opts)
%
% Inputs:
%   ecg   - ECG signal (vector), ideally filtered with positive polarity
%           (the signal output of get_RR)
%   peaks - R-peak sample indices (NaN and beats too close to the edges
%           are ignored)
%   fs    - sampling rate (Hz)
%   opts  - optional struct:
%           .win_ms     - beat epoch (ms, default [-200 400])
%           .amp_win_ms - window for the R amplitude (ms, default [-25 25])
%           .r_thresh   - low-correlation threshold (default 0.7)
%
% Outputs:
%   qa - struct:
%        .NBeats      - number of usable beats
%        .MedCorr     - median template correlation
%        .FracLowCorr - fraction of beats with r < r_thresh   (flag 'shape' if > 0.10)
%        .FracLowAmp  - fraction of low-amplitude beats       (flag 'amp'   if > 0.10)
%        .RRlag1      - RR lag-1 autocorrelation, NaN if <= 10 intervals
%                                                             (flag 'altRR' if < -0.3)
%        .MedRR       - median RR (ms)             (flag 'rate' if < 500 or > 1500)
%        .flag        - 'ok', the flags joined by '+', or 'too-few-beats'
%                       (< 20 usable beats; other fields then NaN)
%
% Ported from the Heartbeam project (HEP_neurofeedback/functions/rpeak_qa.m),
% where it flagged 8 of 78 sessions, including T-wave locking the detector
% itself produced (18.5% low-amplitude beats, RR lag-1 = -0.31).
%
% Copyright (C) - Cedric Cannard, 2026

function qa = rpeak_qa(ecg, peaks, fs, opts)

if nargin < 4 || isempty(opts), opts = struct; end
win_ms     = getdef(opts, 'win_ms',     [-200 400]);
amp_win_ms = getdef(opts, 'amp_win_ms', [-25 25]);
r_thresh   = getdef(opts, 'r_thresh',   0.7);

qa = struct('NBeats',NaN, 'MedCorr',NaN, 'FracLowCorr',NaN, 'FracLowAmp',NaN, ...
    'RRlag1',NaN, 'MedRR',NaN, 'flag','too-few-beats');

sig  = double(ecg(:)');
nPts = numel(sig);
pre  = round(abs(win_ms(1))/1000*fs);
post = round(win_ms(2)/1000*fs);
rp   = round(peaks(:)');
rp   = rp(~isnan(rp) & rp > pre & rp <= nPts - post);
if numel(rp) < 20, return; end

% Beat epochs [nBeats x nTimes]
E  = sig(rp(:) + (-pre:post));
ok = all(isfinite(E), 2);
if sum(ok) < 20, return; end
E  = E(ok,:);
qa.NBeats = size(E,1);

% 1. Template correlation (Pearson r of each beat with the median beat)
tmpl = median(E, 1);
Ec   = E - mean(E, 2);
tc   = tmpl - mean(tmpl);
r    = (Ec * tc(:)) ./ max(sqrt(sum(Ec.^2,2) * sum(tc.^2)), eps);
qa.MedCorr     = median(r);
qa.FracLowCorr = mean(r < r_thresh);

% 2. R amplitude: peak-to-peak within amp_win_ms around t = 0
t_ms  = linspace(win_ms(1), win_ms(2), size(E,2));
amask = t_ms >= amp_win_ms(1) & t_ms <= amp_win_ms(2);
amp   = max(E(:,amask),[],2) - min(E(:,amask),[],2);
qa.FracLowAmp = mean(amp < 0.5*median(amp));

% 3. RR alternation (RR in ms, from all in-bounds peaks)
rr = diff(rp(:))/fs*1000;
qa.MedRR = median(rr);
if numel(rr) > 10
    x = rr - mean(rr);
    qa.RRlag1 = (x(1:end-1)' * x(2:end)) / max(sum(x.^2), eps);
end

% Triage (thresholds from Heartbeam)
f = {};
if qa.FracLowCorr > 0.10,               f{end+1} = 'shape'; end
if qa.FracLowAmp  > 0.10,               f{end+1} = 'amp';   end
if qa.RRlag1      < -0.3,               f{end+1} = 'altRR'; end
if qa.MedRR < 500 || qa.MedRR > 1500,   f{end+1} = 'rate';  end
if isempty(f), qa.flag = 'ok'; else, qa.flag = strjoin(f,'+'); end
end

% Field f of struct s, or default d if missing or empty
function v = getdef(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
