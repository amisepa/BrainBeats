% Returns the starting time (in seconds) of each window to be analyzed.
%
% INPUT:
%   t: a single row of time of the rr interval (in s)
%
% OUTPUT:
%   w : array with starting timeof each window to be analyzed (in s)

function w = get_sqiwindow(t)

increment = 1;
windowlength = 10;

nx = floor(t(end));               % length of sequence
overlap = windowlength-increment;   % number of overlapping elements
Nwinds = fix((nx-overlap)/(windowlength-overlap));    % number of sliding windows
% Initialize output matrix
w = (0:(Nwinds-1))*(windowlength-overlap);  % starting index of each windows

end