% CFA_REPORT - Print and store the heart-locked EEG amplitude before and after
% a heart artifact removal.
%
% Usage:
%   S = cfa_report(before, after, label)
%
% Inputs:
%   before, after - [QRS late] amplitudes from CFA_AMPLITUDE
%   label         - what was done, for the printed line
%
% Output:
%   S - struct: cfa_before, cfa_after (QRS window), cfa_late_before,
%       cfa_late_after (late window)
%
% Copyright (C) - Cedric Cannard, 2026

function S = cfa_report(before, after, label)
S = struct('cfa_before',before(1), 'cfa_after',after(1), ...
    'cfa_late_before',before(2), 'cfa_late_after',after(2));
fprintf('Heart-locked EEG amplitude (%s): -50 to 100 ms %.2f -> %.2f uV (%.0f%% reduction); 200 to 500 ms %.2f -> %.2f uV (%.0f%%). \n', ...
    label, before(1), after(1), 100*(1 - after(1)/before(1)), before(2), after(2), 100*(1 - after(2)/before(2)));
