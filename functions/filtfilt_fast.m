function X = filtfilt_fast(varargin)
% FILTFILT_FAST - Like filtfilt(), but faster when filter and signal are long (and A = 1).
%
% Usage:
%   Y = filtfilt_fast(B, A, X)
%   Y = filtfilt_fast(N, F, A, X)
%
% Inputs:
%   B, A - filter coefficients (A ~= 1 falls back to filtfilt)
%   X    - signal (time x channels), filtered along the first dimension
%   N, F, A - alternatively, design an order-N FIR filter (design_fir) with
%          amplitude response A at normalized frequencies F (0 to 1)
% Outputs:
%   Y    - zero-phase filtered signal (same size and class as X)
%
% Filters forward and backward with filter_fast (FFT convolution for long
% signals and filters), after reflecting the signal at both ends.
%
% See also:
%   filtfilt, filter
% 
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-07-14

% Copyright (C) Christian Kothe, SCCN, 2010, ckothe@ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

if nargin == 3
    [B A X] = deal(varargin{:});
elseif nargin == 4
    [N F M X] = deal(varargin{:});
    B = design_fir(N,F,sqrt(M)); A = 1; % note: we use the sqrt() because we run forward and backward
else
    help filtfilt_fast;
    return;
end

if A == 1
    was_single = strcmp(class(X),'single');
    w = length(B); t = size(X,1);    
    % extrapolate length(B) samples at each end (point reflection) to limit edge effects
    X = double([bsxfun(@minus,2*X(1,:),X(1+mod(((w+1):-1:2)-1,t),:)); X; bsxfun(@minus,2*X(t,:),X(1+mod(((t-1):-1:(t-w))-1,t),:))]);
    % filter, reverse
    X = filter_fast(B,A,X); X = X(length(X):-1:1,:);
    % filter, reverse
    X = filter_fast(B,A,X); X = X(length(X):-1:1,:);
    % remove extrapolated pieces
    X([1:w t+w+(1:w)],:) = [];
    if was_single
        X = single(X); end    
else    
    % fall back to filtfilt for the IIR case
    X = filtfilt(B,A,X);
end
