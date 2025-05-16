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

% function [NN, t_NN, flagged_beats] = clean_rr(t_rr, rr, params)
% 
% %% Parameters
% fs = params.fs;
% if isfield(params,'rr_physlimlow') && ~isempty(params.rr_physlimlow)
%     lowerphysiolim = params.rr_physlimlow;
% else
%     lowerphysiolim = .375;    % default = .375 s (healthy range = .6-1 s)
% end
% if isfield(params,'rr_physlimhigh') && ~isempty(params.rr_physlimhigh)
%     upperphysiolim = params.rr_physlimhigh;  % (healthy range = .6-1 s)
% else
%     upperphysiolim = 2;    % default = 2 s
% end
% if isfield(params,'rr_gaplim') && ~isempty(params.rr_gaplim)
%     gaplimit = params.rr_gaplim;
% else
%     gaplimit = 2;    % default = 2 s 
% end
% if isfield(params,'rr_changelim') && ~isempty(params.rr_changelim)
%     changeLimit = params.rr_changelim;
% else
%     changeLimit = .2;    % default = .2
% end
% if isfield(params,'rr_correct') && ~isempty(params.rr_correct)
%     interpMeth = params.rr_correct;
% else
%     interpMeth = 'pchip';    % default = 'pchip'
% end
% 
% 
% %% Run
% 
% % prep
% % t_rr(1) = [];  % remove 1st heartbeat
% Rpeaks = repmat('N', [length(rr) 1]);
% 
% % Remove data that are too close together (not counted in total signal removed)
% idx_remove = find(diff(t_rr) < 1/fs);
% rr(idx_remove+1) = [];
% t_rr(idx_remove) = [];
% Rpeaks(idx_remove) = [];
% clear idx_remove;
% 
% % Remove artifact Rpeaks (not counted in total signal removed)
% idx_remove = strcmp('|', Rpeaks); % find artifacts
% rr(idx_remove) = [];
% Rpeaks(idx_remove) = [];
% t_rr(idx_remove) = [];
% 
% % find artifacts
% idx_remove2 = strcmp('~', Rpeaks);
% rr(idx_remove2) = [];
% Rpeaks(idx_remove2) = [];
% t_rr(idx_remove2) = [];
% clear idx_remove idx_remove2
% 
% % Remove large RR intervals caused by gaps (not counted in total signal removed)
% idx_remove = find(rr >= gaplimit);
% rr(idx_remove) = [];
% t_rr(idx_remove) = [];
% Rpeaks(idx_remove) = [];
% clear idx_remove;
% 
% % Find non-R peaks
% goodbeats = strcmp(cellstr(Rpeaks),'N');
% badbeats = ~goodbeats;
% outliers = false(length(badbeats),1);
% 
% % remove beats after non-R beats
% for i = 1:length(badbeats)
%     if badbeats(i)
%         outliers(i) = 1;
%         outliers(i+1) = 1;
%     end
% end
% 
% % remove extra points that may have been introduced at the end of the file
% if (length(outliers) > length(badbeats))
%     z = length(outliers);
%     outliers(z) = [];
% end
% if length(Rpeaks) ~= length(rr)
%     outliers = zeros(length(rr),1);
% end
% 
% % Find RR over given percentage change
% idxRRtoBeRemoved = FindSpikesInRR(rr, changeLimit);
% 
% % Combine Rpeaks and percentage outliers
% outliers_combined = idxRRtoBeRemoved(:) + outliers(:);
% outliers = logical(outliers_combined);
% 
% % Remove extra points that may have been introduced at the end of the file
% if length(outliers) > length(outliers)
%     fprintf('Too many Outliers - Removing last Point\n');
%     z = length(outliers);
%     outliers(z) = [];
% end
% 
% % Remove or interpolate outliers
% idx_outliers = find(outliers == 1);
% numOutliers = length(idx_outliers);
% % rr_original = rr;
% rr(idx_outliers) = NaN;
% if strcmp(interpMeth, 'remove')
%     NN_Outliers = rr;
%     NN_Outliers(idx_outliers) = [];
%     t_Outliers = t_rr;
%     t_Outliers(idx_outliers) = [];
% else
%     NN_Outliers = interp1(t_rr,rr,t_rr,interpMeth);
%     t_Outliers = t_rr;
% end
% 
% % Identify non-physiologic beats
% toohigh = NN_Outliers > upperphysiolim;
% toolow = NN_Outliers < lowerphysiolim;
% idx_toolow = find(toolow == 1);
% NN_NonPhysBeats = NN_Outliers;
% NN_NonPhysBeats(idx_toolow) = NaN;
% numOutliers = numOutliers + length(idx_toolow);
% if strcmp(interpMeth, 'remove')
%     NN_NonPhysBeats(idx_toolow) = [];
%     t_NonPhysBeats = t_Outliers;
%     t_NonPhysBeats(idx_toolow) = [];
% else
%     NN_NonPhysBeats = interp1(t_Outliers,NN_NonPhysBeats,t_Outliers,interpMeth);
%     t_NonPhysBeats = t_Outliers;
%     flagged_beats = logical(outliers(:) + toohigh(:) + toolow(:));
% end
% 
% % Interpolate beats that are too Fast
% toohigh = NN_NonPhysBeats > upperphysiolim;
% idx_outliers_2ndPass = find(logical(toohigh(:)) ~= 0);
% NN_TooFastBeats = NN_NonPhysBeats;
% NN_TooFastBeats(idx_outliers_2ndPass) = NaN;
% numOutliers = numOutliers + length(idx_outliers_2ndPass);
% if strcmp(interpMeth, 'remove')
%     flagged_beats = numOutliers;
%     NN_TooFastBeats(idx_outliers_2ndPass) = [];
%     t_TooFastBeats = t_NonPhysBeats;
%     t_TooFastBeats(idx_outliers_2ndPass) = [];
% else
%     NN_TooFastBeats = interp1(t_NonPhysBeats,NN_TooFastBeats,t_NonPhysBeats,'spline','extrap');
%     % NN_TooFastBeats = interp1(t_NonPhysBeats,NN_TooFastBeats,t_NonPhysBeats,interpMeth,'extrap');
%     t_TooFastBeats = t_NonPhysBeats;
% end
% 
% % if params.vis_cleaning
% %     figure('color','w');
% %     plot(t_rr,rr_original,t_Outliers,NN_Outliers);
% %     legend('raw','interp1(after outliers removed)')
% %     hold on
% %     plot(t_NonPhysBeats,NN_NonPhysBeats+.01);
% %     hold on; plot(t_NonPhysBeats,toolow,'o')
% %     legend('raw','interp1(after outliers removed)',...
% %         'interp2(after too low)','toolow')
% %     hold on; plot(t_TooFasyBeats,NN_TooFastBeats);
% %     legend('raw','interp1(after outliers removed)',...
% %         'interp2(after too low)','toolow','interp3 (after too fast removed)')
% % end
% 
% % Remove erroneous data at the end of a record (i.e. a un-physiologic point
% % caused by removing data at the end of a record)
% while NN_TooFastBeats(end) > upperphysiolim
%     NN_TooFastBeats(end) = [];  t_TooFastBeats(end) = [];
% end
% 
% NN = NN_TooFastBeats;
% t_NN = t_TooFastBeats;
% 
% end
% 
% %% clean RR intervals that change more than a given threshold
% % (eg., th = 0.2 = 20%) with respect to the median value of the previous 5
% % and next 5 RR intervals (using a forward-backward approach).
% %
% % INPUTS:
% %       RR : a single row of rr interval data in seconds
% %       th : threshold percent limit of change from one interval to the next
% % OUTPUTS:
% %       idxRRtoBeRemoved : a single vector of indexes related to RR
% %                          intervals corresponding to a change > th
% 
% function idxRRtoBeRemoved = FindSpikesInRR(RR, th)
% 
% if size(RR,1)>size(RR,2)
%     RR = RR';
% end
% 
% % Forward search
% FiveRR_MedianVal = medfilt1(RR,5); % compute as median RR(-i-2: i+2)
% 
% % shift of three position to align with to corresponding RR
% FiveRR_MedianVal = [RR(1:5) FiveRR_MedianVal(3:end-3)];
% rr_above_th = (abs(RR-FiveRR_MedianVal)./FiveRR_MedianVal)>=th;
% 
% RR_forward = RR;
% RR_forward(rr_above_th) = NaN;
% 
% % Backward search
% RRfilpped = fliplr(RR);
% FiveRR_MedianVal = medfilt1(RRfilpped,5); % compute as median RR(-i-2: i+2)
% % shift of three position to aligne with to corresponding RR
% FiveRR_MedianVal = [RRfilpped(1:5) FiveRR_MedianVal(3:end-3)];
% rr_above_th = find(abs(RRfilpped-FiveRR_MedianVal)./FiveRR_MedianVal>=th);
% rr_above_th = sort(length(RR)-rr_above_th+1);
% 
% RR_backward = RRfilpped;
% RR_backward(rr_above_th) = NaN;
% 
% % Combine
% idxRRtoBeRemoved = (isnan(RR_forward) & isnan(RR_backward));
% 
% end



function [nn_intervals, nn_t, idx_corrected] = clean_rr(rr_t, rr_intervals)
    % clean_rr_intervals - Detects and corrects abnormal heartbeats in RR interval time series.
    %
    % INPUT:
    % rr_t         - Vector of timestamps (seconds) corresponding to each RR interval.
    % rr_intervals - Vector of RR intervals in seconds.
    % peaks        - Vector of sample indices of detected peaks.
    %
    % OUTPUT:
    % nn_intervals - Vector of corrected NN intervals in seconds.
    % new_peaks    - Updated sample indices of peaks after corrections.
    % idx_corrected - Logical array indicating which intervals were corrected (true if corrected).
    
    % Parameters
    lower_limit = 0.375;  % Minimum valid RR interval (0.375 s = 160 bpm)
    upper_limit = 1.5;      % Maximum valid RR interval (1.5 s = 40 bpm and 0.66 hz)
    
    nn_intervals = rr_intervals;
    nn_t = rr_t';
    % new_peaks = peaks;
    idx_corrected = false(length(rr_intervals), 1);

    % Step 1: Remove overly short RR intervals
    % short_rr = diff(rr_t) <= lower_limit;
    short_rr = rr_intervals <= lower_limit;
    if any(short_rr)
        warning('Removing %g overly short RR intervals', sum(short_rr));
    end
    nn_intervals(short_rr) = [];
    nn_t(short_rr) = [];
    % new_peaks(short_rr) = [];
    idx_corrected(short_rr) = true;
    
    % figure('Color', 'w'); hold on;
    % plot(nn_t, nn_intervals, 'Color',[0.6350 0.0780 0.1840], 'LineWidth', 2, 'DisplayName', 'Short RR');

    % Step 2: Recalculate gaps based on updated nn_intervals
    gaps = nn_intervals > upper_limit;
    if any(gaps)
        warning('Interpolating %g gaps in the RR intervals', sum(gaps));
    end
    idx_corrected(gaps) = true;
    nn_intervals = interp1(nn_t(~gaps), nn_intervals(~gaps), nn_t, 'spline', 'extrap');
    % plot(nn_t, nn_intervals, 'Color',[0.3010 0.7450 0.9330], 'LineWidth', 2, 'DisplayName', 'Gaps');

    % Step 3: Interpolate outliers
    outliers = isoutlier(nn_intervals, 'median');
    if any(outliers)
        warning('Interpolating %g outlier RR intervals', sum(outliers));
    end
    idx_corrected(outliers) = true;
    % plot(nn_t, nn_intervals, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2, 'DisplayName', 'Outliers');
    % nn_intervals = interp1(nn_t(~outliers), nn_intervals(~outliers), nn_t, 'pchip');
    nn_intervals = interp1(nn_t(~outliers), nn_intervals(~outliers), nn_t, 'spline', 'extrap');

    % Step 4L interpolate sharp spikes
    spikes = FindSpikesInRR(nn_intervals, .2);
    if any(spikes)
        warning('Interpolating %g outlier RR intervals', sum(outliers));
    end
    idx_corrected(spikes) = true;
    % plot(nn_t, nn_intervals, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2, 'DisplayName', 'Outliers');
    % nn_intervals = interp1(nn_t(~spikes), nn_intervals(~spikes), nn_t, 'pchip');
    nn_intervals = interp1(nn_t(~spikes), nn_intervals(~spikes), nn_t, 'spline', 'extrap');

    % plot(nn_t, nn_intervals, 'k', 'LineWidth', 2, 'DisplayName', 'Final');
    % legend('show');

    % Output results
    if sum(idx_corrected)==0
        idx_corrected = [];
        disp("No abnormal RR intervals detected.")
    else
        fprintf('Abnromal RR intervals detected and corrected: %g\n', sum(idx_corrected));
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