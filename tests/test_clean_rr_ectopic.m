function test_clean_rr_ectopic
% TEST_CLEAN_RR_ECTOPIC - clean_rr removes ectopic beats (premature beat at
% 0.7 x the sinus interval, full compensatory pause) and leaves sinus rhythm
% with respiratory sinus arrhythmia untouched, with the 'adaptive' rule
% (default since v1.7) and the 'ratio' rule. No EEGLAB needed.
%
% Both rules on real recordings: tests/compare_ectopic_rules.m.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'functions'));
rng(1); fs = 250;
t = 0; bt = 0;
while t < 300                                   % sinus: 800 ms, RSA 40 ms at 0.25 Hz
    t = t + 0.8 + 0.04*sin(2*pi*0.25*t) + 0.01*randn;
    bt(end+1,1) = t; %#ok<AGROW>
end
rmssd = @(x) sqrt(mean(diff(x*1000).^2));
idx = round(linspace(30, numel(bt)-30, 6));
bt2 = bt;
for i = idx
    bt2(i+1) = bt(i) + 0.7*(bt(i+1)-bt(i));     % ectopic replaces sinus beat i+1
end
for rule = {'adaptive' 'ratio'}
    [nn0, ~, ~, ~, ii0] = run(bt, fs, rule{1});
    assert(numel(nn0) == numel(bt)-1 && ~any(ii0), '%s: clean_rr changed a sinus series', rule{1})
    [nnE, ~, ~, ~, iiE] = run(bt2, fs, 'none');
    [nn, ~, ~, ~, ii] = run(bt2, fs, rule{1});
    assert(sum(ii) == 6, '%s: %d of 6 ectopic beats replaced', rule{1}, sum(ii))
    assert(~any(iiE), 'without the ectopic rule, beats were replaced')
    r0 = rmssd(nn0); rE = rmssd(nnE); r = rmssd(nn);
    assert(abs(r/r0 - 1) < 0.1, '%s: RMSSD %.1f ms after cleaning vs %.1f ms without ectopic beats', rule{1}, r, r0)
    fprintf('clean_rr ectopic OK (%s): RMSSD %.1f ms (sinus), %.1f with 6 ectopic beats, %.1f after cleaning\n', rule{1}, r0, rE, r);
end

% RR_ALTERNATION: none on sinus series (strong respiratory sinus arrhythmia
% included) or once occasional ectopic beats are replaced; above the 10%
% warning of BRAINBEATS_PROCESS in bigeminy, which the cleaning only partly
% corrects
nnR = 0.8 + 0.15*sin(2*pi*(1:300)'*0.8/4);       % +/-19%, 4-s breaths
bg = repmat([0.6; 1.0], 150, 1) + 0.01*randn(300, 1);   % bigeminy: premature beat, pause
[nnB, ~, ~, ~, iiB] = run(cumsum([0; bg]), fs, 'adaptive');
a = [rr_alternation(nn0) rr_alternation(nn) rr_alternation(nnR) rr_alternation(nnB)];
assert(all(a(1:3) == 0), 'rr_alternation on sinus series: %s', mat2str(a(1:3)))
assert(a(4) > 0.1, 'rr_alternation in bigeminy: %.2f (%d beats replaced)', a(4), sum(iiB))
fprintf('rr_alternation OK: 0 on sinus and cleaned series, %.2f in bigeminy\n', a(4));
end

function [nn, nn_t, pk, bad, interp] = run(bt, fs, rule)
RR = diff(bt); RR_t = bt(2:end); pk = round(RR_t*fs) + 1;
[~, nn, nn_t, pk, bad, interp] = evalc('clean_rr(RR_t, RR, [], pk, ''ectopic_rule'', rule)');
end
