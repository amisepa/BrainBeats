function test_clean_rr_ectopic
% TEST_CLEAN_RR_ECTOPIC - clean_rr removes ectopic beats (premature beat at
% 0.7 x the sinus interval, full compensatory pause) and leaves sinus rhythm
% with respiratory sinus arrhythmia untouched. No EEGLAB needed.
%
% On 6 ds003838 recordings (2026-11-03, rule at 20%): no beat flagged in the 5
% clean ones; 5 in sub-038, all abnormal on the ECG. Injected premature beats
% at 0.70 x were detected in 1-6 of 6 per recording, at 0.75 x in 1-4: the
% rule is conservative where the sinus rhythm varies a lot.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'functions'));
rng(1); fs = 250;
t = 0; bt = 0;
while t < 300                                   % sinus: 800 ms, RSA 40 ms at 0.25 Hz
    t = t + 0.8 + 0.04*sin(2*pi*0.25*t) + 0.01*randn;
    bt(end+1,1) = t; %#ok<AGROW>
end
rmssd = @(x) sqrt(mean(diff(x*1000).^2));
[nn0, ~, ~, ~, ii0] = run(bt, fs);
assert(numel(nn0) == numel(bt)-1 && ~any(ii0), 'clean_rr changed a sinus series')

idx = round(linspace(30, numel(bt)-30, 6));
bt2 = bt;
for i = idx
    bt2(i+1) = bt(i) + 0.7*(bt(i+1)-bt(i));     % ectopic replaces sinus beat i+1
end
[nnE, ~, ~, ~, iiE] = run(bt2, fs, Inf);
[nn, ~, ~, ~, ii] = run(bt2, fs);
assert(sum(ii) == 6, '%d of 6 ectopic beats replaced', sum(ii))
assert(~any(iiE), 'without the ectopic rule, beats were replaced')
r0 = rmssd(nn0); rE = rmssd(nnE); r = rmssd(nn);
assert(abs(r/r0 - 1) < 0.1, 'RMSSD %.1f ms after cleaning vs %.1f ms without ectopic beats', r, r0)
fprintf('clean_rr ectopic OK: RMSSD %.1f ms (sinus), %.1f with 6 ectopic beats, %.1f after cleaning\n', r0, rE, r);
end

function [nn, nn_t, pk, bad, interp] = run(bt, fs, ratio)
RR = diff(bt); RR_t = bt(2:end); pk = round(RR_t*fs) + 1;
if nargin < 3, ratio = 0.2; end
[~, nn, nn_t, pk, bad, interp] = evalc('clean_rr(RR_t, RR, [], pk, ''ectopic_ratio'', ratio)');
end
