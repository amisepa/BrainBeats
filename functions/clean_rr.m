%% Preprocesses RR interval data for running HRV analyses.
% Noise and non-normal beats are either removed from the RR interval data
% or interpolated using a interpolation method of choice.
% The default extrapolation behavior is 'extrap' for 'spline', 'pchip' and 'makima',
% and EXTRAPVAL = NaN (NaN+NaN*1i for complex values) for the other methods.

% INPUTS:
%   t_rr - time of the rr interval data
%   rr - a single row of rr interval data (in s)
%   vis - plot (1) or not (0)
%
% OUTPUTS:
%   nn - normal normal interval data
%   t_nn - time of NN interval
%   fs - sample rate (Hz)
%   flagged_beats - percent of original data removed
%
% PLEASE CITE:
%   Clifford, G. (2002). "Characterizing Artefact in the Normal
%   Human 24-Hour RR Time Series to Aid Identification and Artificial
%   Replication of Circadian Variations in Human Beat to Beat Heart
%   Rate using a Simple Threshold."
%
%   Vest et al. (2018) "An Open Source Benchmarked HRV Toolbox for Cardiovascular
%   Waveform and Interval Analysis" Physiological Measurement.
%
%   ORIGINAL SOURCE AND AUTHORS:
%       Adriana N. Vest and various authors (Physionet Cardiovascular Signal toolbox).
%
% Copyright (C), BrainBeats, Cedric Cannard, 2022


function [nn_intervals, nn_t, idx_rem, idx_corr] = clean_rr(rr_t, rr_intervals, peak_amp)
    % clean_rr_intervals - Detects and corrects abnormal heartbeats in RR interval time series.
    %
    % INPUT:
    %   rr_t          - Vector of timestamps (seconds) corresponding to each RR interval.
    %   rr_intervals  - Vector of RR intervals in seconds.
    %   peak_amp      - Vector of signal amplitude at detected peaks.
    %
    % OUTPUT:
    %   nn_intervals  - Vector of corrected NN intervals in seconds.
    %   new_peaks     - Updated sample indices of peaks after corrections.
    %   idx_removed   - Logical array indicating which intervals were removed (true if corrected).
    %   idx_corrected - Logical array indicating which intervals were corrected (true if corrected).

    % Parameters
    lower_limit = 0.375;  % Minimum valid RR interval (0.375 s = 160 bpm)
    upper_limit = 1.5;      % Maximum valid RR interval (1.5 s = 40 bpm and 0.66 hz)
    
    % initialize variables
    nn_intervals = rr_intervals;
    nn_t = rr_t';
    idx_rem = false(length(rr_intervals), 1);

    % Step 1: Remove overly short RR intervals
    short_rr = rr_intervals <= lower_limit;
    if any(short_rr)
        warning('Removing %g overly short RR intervals', sum(short_rr));
    end
    nn_intervals(short_rr) = [];
    nn_t(short_rr) = [];
    % new_peaks(short_rr) = [];
    peak_amp(short_rr) = [];
    idx_rem(short_rr) = true;

    % figure('Color', 'w'); hold on;
    % plot(nn_t, nn_intervals, 'Color',[0.6350 0.0780 0.1840], 'LineWidth', 2, 'DisplayName', 'Short RR');

    % Step 2: Recalculate gaps based on updated nn_intervals
    idx_corr = false(length(nn_intervals), 1);
    gaps = nn_intervals > upper_limit;
    if any(gaps)
        warning('Interpolating %g gaps in the RR intervals', sum(gaps));
    end
    idx_corr(gaps) = true;
    nn_intervals = interp1(nn_t(~gaps), nn_intervals(~gaps), nn_t, 'spline', 'extrap');
    % plot(nn_t, nn_intervals, 'Color',[0.3010 0.7450 0.9330], 'LineWidth', 2, 'DisplayName', 'Gaps');

    % Step 3: Interpolate outliers
    outliers = isoutlier(nn_intervals, 'mean') | isoutlier(peak_amp,'mean');
    if any(outliers)
        warning('Interpolating %g outlier RR intervals', sum(outliers));
    end
    idx_corr(outliers) = true;
    % plot(nn_t, nn_intervals, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2, 'DisplayName', 'Outliers');
    % nn_intervals = interp1(nn_t(~outliers), nn_intervals(~outliers), nn_t, 'pchip');
    nn_intervals = interp1(nn_t(~outliers), nn_intervals(~outliers), nn_t, 'spline', 'extrap');

    % Step 4: interpolate sharp spikes
    spikes = FindSpikesInRR(nn_intervals, .2);
    if any(spikes)
        warning('Interpolating %g outlier RR intervals', sum(outliers));
    end
    idx_corr(spikes) = true;
    % plot(nn_t, nn_intervals, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2, 'DisplayName', 'Outliers');
    % nn_intervals = interp1(nn_t(~spikes), nn_intervals(~spikes), nn_t, 'pchip');
    nn_intervals = interp1(nn_t(~spikes), nn_intervals(~spikes), nn_t, 'spline', 'extrap');

    % plot(nn_t, nn_intervals, 'k', 'LineWidth', 2, 'DisplayName', 'Final');
    % legend('show');

    % Output results
    if sum(idx_corr)==0
        idx_corr = [];
        disp("No abnormal RR intervals detected.")
    else
        fprintf('Abnromal RR intervals detected and corrected: %g\n', sum(idx_corr));
    end


end

%% clean RR intervals that change more than a given threshold
% (eg., th = 0.2 = 20%) with respect to the median value of the previous 5
% and next 5 RR intervals (using a forward-backward approach).
%
% INPUTS:
%       RR : a single row of rr interval data in seconds
%       th : threshold percent limit of change from one interval to the next
% OUTPUTS:
%       idxRRtoBeRemoved : a single vector of indexes related to RR
%                          intervals corresponding to a change > th

function idxRRtoBeRemoved = FindSpikesInRR(RR, th)

if size(RR,1)>size(RR,2)
    RR = RR';
end

% Forward search
FiveRR_MedianVal = medfilt1(RR,5); % compute as median RR(-i-2: i+2)

% shift of three position to align with to corresponding RR
FiveRR_MedianVal = [RR(1:5) FiveRR_MedianVal(3:end-3)];
rr_above_th = (abs(RR-FiveRR_MedianVal)./FiveRR_MedianVal)>=th;

RR_forward = RR;
RR_forward(rr_above_th) = NaN;

% Backward search
RRfilpped = fliplr(RR);
FiveRR_MedianVal = medfilt1(RRfilpped,5); % compute as median RR(-i-2: i+2)
% shift of three position to aligne with to corresponding RR
FiveRR_MedianVal = [RRfilpped(1:5) FiveRR_MedianVal(3:end-3)];
rr_above_th = find(abs(RRfilpped-FiveRR_MedianVal)./FiveRR_MedianVal>=th);
rr_above_th = sort(length(RR)-rr_above_th+1);

RR_backward = RRfilpped;
RR_backward(rr_above_th) = NaN;

% Combine
idxRRtoBeRemoved = (isnan(RR_forward) & isnan(RR_backward));

end