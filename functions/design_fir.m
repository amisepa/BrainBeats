function B = design_fir(N,F,A,nfft,W)
% DESIGN_FIR - Design an FIR filter with the frequency-sampling method.
%
% The amplitude response is interpolated (piecewise cubic, pchip) between the
% specified frequency points, given linear phase, and windowed.
%
% Usage:
%   B = design_fir(N, F, A, nFFT, W)
%
% Inputs:
%   N    - filter order
%   F    - frequencies at which amplitudes are defined, normalized to Nyquist
%          (from 0 to 1; avoid too sharp transitions)
%   A    - amplitudes, one per frequency in F
%   nFFT - (optional) number of FFT bins (default max(512, next power of 2 >= N))
%   W    - (optional) window, N+1 samples (default Hamming)
% Outputs:
%   B    - filter kernel (1 x N+1)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-08-14

% Copyright (C) Christian Kothe, SCCN, 2013, ckothe@ucsd.edu
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

if nargin < 4 || isempty(nfft)
    nfft = max(512,2^ceil(log(N)/log(2))); end
if nargin < 5
    W = 0.54 - 0.46*cos(2*pi*(0:N)/N); end

% calculate interpolated frequency response
F = interp1(round(F*nfft),A,(0:nfft),'pchip');

% set phase & transform into time domain
F = F .* exp(-(0.5*N)*sqrt(-1)*pi*(0:nfft)./nfft);
B = real(ifft([F conj(F(end-1:-1:2))]));

% apply window to kernel
B = B(1:N+1).*W(:)';
