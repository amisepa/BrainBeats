% BB_FAIL - Report a failed input check: an error on the command line, so that
% scripts stop, or a dialog in the GUI (the caller then returns).
%
% Usage:
%   bb_fail(isCmd, msg)
%
% Copyright (C) - Cedric Cannard, 2026

function bb_fail(isCmd, msg)
if isCmd
    error('BrainBeats:input', '%s', char(msg));
end
errordlg(char(msg), 'BrainBeats');
