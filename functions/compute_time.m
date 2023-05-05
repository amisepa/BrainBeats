function time = compute_time(NN, params)

% if params.winSize < 300
%     warning('File length is too short to reliably estimate HRV time-domain metrics. At least 5 minutes of data are recommended (see Shaffer and Ginsberg, 2017).')
% end

% NN metrics
if params.hrv_time
    time.NN_mean = mean(NN.*1000);      % in ms
    time.NN_median = median(NN.*1000);  % in ms
    time.NN_mode = mode(NN.*1000);      % in ms
    time.NN_var = var(NN.*1000);        % in ms
    time.NN_skew = skewness(NN);
    time.NN_kurt = kurtosis(NN);
    time.NN_iqr = iqr(NN.*1000);        % in ms

    % SDNN
    time.SDNN = std(NN.*1000);   % in ms

    % RMSSD (sqrt of the mean squared time diff between heartbeats)
    time.RMSSD = sqrt(mean(diff(NN.*1000).^2));  % in ms

    % pNN50 (fraction of differences larger than alpha = 50)
    alpha = 50;
    time.pNN50 = sum( abs(diff(NN)) >= alpha/1000 )/length(diff(NN));

end

