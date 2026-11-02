% BB_WARNDLG - Warning dialog, except under 'matlab -batch', where nobody can
% close it (the message is then only in the warning printed by the caller).
%
% Usage:
%   bb_warndlg(msg)
%   bb_warndlg(msg, title)
%
% Copyright (C) - Cedric Cannard, 2026

function bb_warndlg(msg, title)
if exist('batchStartupOptionUsed') > 0 && batchStartupOptionUsed %#ok<EXIST>
    return
end
if nargin < 2, title = 'BrainBeats'; end
warndlg(msg, title);
