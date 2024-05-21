%% Resample NN intervals to compute PSD with pwelch or FFT
%
%   INPUT:
%       win_idx         - start and end time of the NN interval
%       NN              - vector of NN intervals to be resampled
%       sf              - 7
%       interp_method   - 'cub' or 'spline' recommended
%
%   OUTPUT
%       NN_resamp      - resampled NN intervals

function [NN_resamp, t_resamp] = resample_NN(NN_times,NN,sf,interp_method)

% time index
t_resamp = NN_times(1):1/sf:NN_times(end);

% Resample with interpolation method of choice
switch interp_method
    case 'cub'
        NN_resamp = interp1(NN_times,NN,t_resamp','spline')'; % cubic spline interpolation (default)
    case 'lin'
        NN_resamp = interp1(NN_times,NN,t_resamp','linear')'; % linear interpolation
end
