% RESAMPLE_NN - Resample an NN interval series on a regular time grid (e.g.
% for PSD estimation with pwelch or FFT).
%
% Usage:
%   [NN_resamp, t_resamp] = resample_NN(NN_times, NN, sf, interp_method)
%
% Inputs:
%   NN_times      - time of each NN interval (s)
%   NN            - NN intervals (s), same length as NN_times
%   sf            - resampling frequency (Hz; e.g. 7 Hz for HRV spectra)
%   interp_method - 'cub' or 'spline' (cubic spline), 'lin' or 'linear',
%                   or 'pchip' (shape-preserving cubic); no default
%
% Outputs:
%   NN_resamp - resampled NN intervals (s), row
%   t_resamp  - resampled time grid (s), row, from NN_times(1) to NN_times(end)
%
% Copyright (C) - Cedric Cannard, 2024

function [NN_resamp, t_resamp] = resample_NN(NN_times,NN,sf,interp_method)

% regular time grid
t_resamp = NN_times(1):1/sf:NN_times(end);

% Resample with interpolation method of choice
switch lower(interp_method)
    case {'cub' 'spline'}
        NN_resamp = interp1(NN_times,NN,t_resamp','spline')'; % cubic spline interpolation
    case {'lin' 'linear'}
        NN_resamp = interp1(NN_times,NN,t_resamp','linear')'; % linear interpolation
    case 'pchip'
        NN_resamp = interp1(NN_times,NN,t_resamp','pchip')';  % shape-preserving cubic
    otherwise
        error("resample_NN: unknown interp_method '%s'. Use 'cub'/'spline', 'lin'/'linear' or 'pchip'.", interp_method)
end
