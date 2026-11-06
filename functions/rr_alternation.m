% RR_ALTERNATION - Share of heartbeats whose interval change reverses at the
% next heartbeat, both changes larger than 15% of the median of the 11
% intervals around them: the short-long alternation of frequent ectopic
% beats and of atrial fibrillation, which CLEAN_RR does not correct.
%
% Usage:
%   f = rr_alternation(nn)
%
% Input:
%   nn - NN intervals (s), e.g. from CLEAN_RR
%
% Output:
%   f  - share of heartbeats (0-1); NaN with fewer than 3 intervals
%
% Notes:
%   Respiratory sinus arrhythmia changes the intervals over several beats
%   and rarely reverses a large change at the next beat. After CLEAN_RR
%   (default rule) on the beats GET_RR detects (2026-11-06):
%   MIT-BIH Arrhythmia Database (tests/validate_ectopic_mitbih.m): 11-44% on
%   the 15 records with frequent ectopic beats, atrial fibrillation or
%   flutter, at most 4.5% on the other 29; 66 resting recordings of healthy
%   adults (ds003838 and the sample data, tests/check_rr_alternation.m): at
%   most 6.9% (sub-055, RMSSD 134 ms), 95th percentile 1.0%. BRAINBEATS_PROCESS
%   warns above 10%.
%
% Copyright (C) - Cedric Cannard, 2026

function f = rr_alternation(nn)
nn = nn(:);
if numel(nn) < 3, f = NaN; return; end
m = movmedian(nn, 11);
d = diff(nn) ./ m(2:end);
f = mean(sign(d(1:end-1)) ~= sign(d(2:end)) & min(abs(d(1:end-1)), abs(d(2:end))) > 0.15);
end
