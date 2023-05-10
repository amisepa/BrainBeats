%% Extract EEG features in time, fequency, and nonlinear domains. 

function EEG = get_eeg_features(signals,times,fs)


% Time domain
EEG.time.trimmed_mean = trimmean(signals,20,2); 
HRV.time.NN_median = median(NN.*1000);  % in ms
HRV.time.NN_mode = mode(NN.*1000);      % in ms
HRV.time.NN_var = var(NN.*1000);        % in ms
HRV.time.NN_skew = skewness(NN);
HRV.time.NN_kurt = kurtosis(NN);
HRV.time.NN_iqr = iqr(NN.*1000);        % in ms







