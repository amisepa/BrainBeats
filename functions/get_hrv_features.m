%% Compute HRV measures on NN time series. 
% 1) Time domain: NN stats, SDNN, RMSSD, pNN50
% 
% 2) Frequency domain: ULF, VLF, LF, HF, LF/HF ratio, TTLPWR.
%   Options: PSD can be calculated using the normalized Lomb-Scargle 
%   periodogram (default), standard Lomb-Scargle periodogram, 
%   Welch, or FFT. Lomb-Scargle periodogram does not require 
%   interpolation or resampling of the data (contrary to welch or FFT), thus 
%   preserving the original information. The Lomb-Scargle method is
%   recommended as it better deals with non-uniformly sampled data, missing 
%   data, noise (common with NN intervals), and does not require resampling. 
%   The normalized version is selected by default (although users can choose 
%   the standard version) by scaling the power by the variance of the signal, 
%   making results more comparable across different recordings or subjects. 
%   If users set hrv_norm to ON, a 2nd level normalization is applied
%   by dividing each band-power by the total power, to provide a more 
%   intuitive measure of the relative contribution of each frequency 
%   component to the overall power. Preferable to do this option when VLF
%   and ULF are available (long time-series required). 
%
% 3) Nonlinear domain: Poincare, fuzzy entropy, fractal dimension, 
%   phase-rectified signal amplitude (PRSA).
%
% Following recommendations by the Task Force of the European Society of 
% Cardiology and the North American Society of Pacing and Electrophysiology
% (1996), minimum data length for each band is 5-10 cycles.
%   - ULF (0-0.003 Hz): at least 24 hours
%   - VLF (0.003-0.04 Hz): at least 28 minutes
%   - LF (0.04-0.15 Hz): at least 125 s
%   - HF (0.15-0.40 Hz): at least 34 s
% 
% To maximize trade-off between time resolution and frequency resolution,
% sliding time windows using these minimum lengths are used for each band.
% Warnings are printed when length is smaller than minimum recommended. 
% 
% Time and nonlinear HRV features are computed over the whole time-series 
% for faster computation and higher reliability. 
% 
% Copyright (C) - Cedric Cannard, 2023

function [HRV, params] = get_hrv_features(NN, NN_times, params)


%% Time domain

if params.hrv_time

    disp('Extracting HRV features in the time domain...')

    % if NN_times(end) < 300
    %     warning('File length is too short to reliably estimate HRV time-domain metrics. At least 5 minutes of data are recommended (see Shaffer and Ginsberg, 2017).')
    % end

    % NN statistics
    % HRV.time.NN_mean = mean(NN.*1000);      % in ms
    % HRV.time.NN_median = median(NN.*1000);  % in ms
    % HRV.time.NN_mode = mode(NN.*1000);      % in ms
    % HRV.time.NN_var = var(NN.*1000);        % in ms
    % HRV.time.NN_skew = skewness(NN);
    % HRV.time.NN_kurt = kurtosis(NN);
    % HRV.time.NN_iqr = iqr(NN.*1000);        % in ms

    % SDNN (standard deviation of the NN intervals)
    HRV.time.SDNN = round(std(NN.*1000),1);   % in ms

    % RMSSD (sqrt of the mean squared time diff between heartbeats)
    HRV.time.RMSSD = round(sqrt(mean(diff(NN.*1000).^2)),1);  % in ms

    % pNN50 (fraction of differences larger than alpha = 50)
    % requires at least 2 min of data (see Ginsberg and Schaffer 2017)
    alpha = 50;
    HRV.time.pNN50 = round(sum( abs(diff(NN)) >= alpha/1000 )/length(diff(NN)),1);

end

%% Frequency domain

if params.hrv_frequency

    disp('Extracting HRV features in the frequency domain...')

    % Parameters
    if isfield(params,'hrv_spec') && ~isempty(params.hrv_spec)
        hrv_spec = params.hrv_spec;
    else
        hrv_spec = 'LombScargle_norm';
        params.hrv_spec_method = 'LombScargle_norm';  % for exportation for users
    end
    if isfield(params,'hrv_norm') && ~isempty(params.hrv_norm)
        norm = params.hrv_norm;
    else
        norm = false;
        params.hrv_normalization = norm;     % for exportation for users
    end
    if isfield(params,'hrv_overlap') && ~isempty(params.hrv_overlap)
        overlap = params.hrv_overlap;
    else
        overlap = .25;  % window overlap (default = 25 %)
        params.hrv_window_overlap = overlap;  % for exportation for users
    end

    % HRV frequency bands (ULF; VLF; LF; HF)
    bands = [ 0 .003; 0.003 .04; .04 .15; 0.15 0.40 ];
    bandNames = {'ULF' 'VLF' 'LF' 'HF'};
    params.hrv_band_freqs = bands;  % for exportation for users
    params.hrv_band_names = bandNames;  % for exportation for users

    % Minimum data length requirements for each band
    % minULF = 86400;       % 24 hours
    % minVLF = 5/0.003;     % 5 cycles/0.03 hz  (in s)
    % minLF = 5/0.04;       % 5 cycles/0.04 hz  (in s)
    % minHF = 5/0.15;       % 5 cycles/0.15 hz  (in s)
    minLength = ceil([ 86400  5/0.003 5/0.04 5/0.15 ]);
    
    for iBand = 1:size(bands,1)

        if NN_times(end) >= minLength(iBand)
    
            % Determine best sliding window length and indices
            winLength = minLength(iBand);
            stepSize = winLength * (1 - overlap);
            nWindows = floor((NN_times(end) - winLength) / stepSize) + 1;
            
            fprintf('Frequency band: %s \n', bandNames{iBand})

            % Compute PSD on each sliding window
            for iWin = 1:nWindows
                fprintf(' - window %g \n', iWin)

                start_idx = (iWin - 1) * stepSize + 1;
                end_idx = start_idx + winLength - 1;
                win_idx = NN_times >= start_idx & NN_times <= end_idx;

                % Frequency resolution and vector for this window
                % nfft = 2^nextpow2(length(NN(win_idx)));    % dynamic nfft based on window length
                nfft = max(2^nextpow2(length(NN(win_idx))), 512);  % use 512 as minimum value for smaller windows
                fvec = bands(iBand,1):1/nfft:bands(iBand,2);
                                
                % Lomb-Scargle Periodogram (no resampling required and best method)
                if strcmp(hrv_spec, 'LombScargle_norm')
                        [pwr,freqs] = plomb(NN(win_idx),NN_times(win_idx),fvec,'normalized'); 
                        fprintf('Computing normalized Lomb-Scargle periodogram on the NN series... \n')
                elseif strcmp(hrv_spec, 'LombScargle')
                        [pwr,freqs] = plomb(NN(win_idx),NN_times(win_idx),fvec,'psd'); 
                        fprintf('Computing standard Lomb-Scargle periodogram on the NN series... \n')
                
                % Welch or FFT (require resampling)
                else
                    % Resample
                    resamp_freq = 7;
                    NN_resamp = resample_NN(NN_times(win_idx),NN(win_idx),resamp_freq,'cub'); % resample RR 
            
                    % Pwelch
                    if strcmp(hrv_spec, 'welch')
                        [pwr,freqs] = pwelch(NN_resamp,minLength(iBand),[],[],resamp_freq);
                        fprintf('Computing pwelch on the NN series... \n')

                    % FFT
                    elseif strcmp(hrv_spec, 'fft')
                        pwr = fft(NN_resamp).*conj(fft(NN_resamp))/length(NN_resamp);
                        freqs = resamp_freq*(0:floor(length(NN_resamp)/2)+1)/length(NN_resamp);
                        pwr = pwr(1:length(freqs));
                        fprintf('Computing FFT on the NN series... \n')
                    end
                end 
                
                % Freq index
                freq_idx = bands(iBand,1) <= freqs & freqs <= bands(iBand,2);
                freq_res = freqs(2)-freqs(1); % resolution

                % Power for each band in ms^2
                if iBand == 1
                    HRV.frequency.ulf(iWin,:) = sum(pwr(freq_idx)*freq_res) * 1e6;      % ULF
                elseif iBand == 2
                    HRV.frequency.vlf(iWin,:) = sum(pwr(freq_idx)*freq_res) * 1e6;      % VLF
                elseif iBand == 3
                    HRV.frequency.lf(iWin,:) = sum(pwr(freq_idx)*freq_res) * 1e6;       % LF
                elseif iBand == 4
                    HRV.frequency.hf(iWin,:) = sum(pwr(freq_idx)*freq_res) * 1e6;       % HF
                end
                
                % Overal power spectra
                PWR{iBand,:} = pwr;         % for exporting overall power spectra
                PWR_freqs{iBand,:} = freqs;  % for exporting overall power spectra

            end
        else
            warning('File length is too short for estimating %s power reliably. At least %1.1f minutes are required. Cannot export this variable.', bandNames{iBand},minLength(iBand)/60)
        end
    end

    % LF/HF ratio (average power across windows)
    if isfield(HRV.frequency,'lf') && isfield(HRV.frequency,'hf')
        HRV.frequency.lfhf = round(mean(HRV.frequency.lf) / mean(HRV.frequency.hf) * 100)/100;
    end
    
    % Total power and normalize if all bands are present (average across windows)
    % if isfield(HRV.frequency,'ulf') && isfield(HRV.frequency,'vlf') ...
    %         && isfield(HRV.frequency,'lf') && isfield(HRV.frequency,'hf')
    if norm
        try
            HRV.frequency.ttlpwr = sum([mean(HRV.frequency.ulf) mean(HRV.frequency.vlf) ...
                    mean(HRV.frequency.lf) mean(HRV.frequency.hf)]);
        catch
            warndlg('HRV total power does not include ULF power (likely due to the short length of the cardiovascular time series).')
            warning('HRV total power does not include ULF power (likely due to the short length of the cardiovascular time series).')
            try
                HRV.frequency.ttlpwr = sum([mean(HRV.frequency.vlf) mean(HRV.frequency.lf) ...
                    mean(HRV.frequency.hf)]);
            catch
                warndlg('HRV total power does not include VLF power (likely due to the short length of the cardiovascular time series).')
                warning('HRV total power does not include VLF power (likely due to the short length of the cardiovascular time series).')
                try
                    HRV.frequency.ttlpwr = sum([mean(HRV.frequency.lf) mean(HRV.frequency.hf)]);
                catch
                    warndlg('Sorry, LF-HRV and HF-HRV power could not be normalized to total power.')
                    warning('Sorry, LF-HRV and HF-HRV power could not be normalized to total power.')
                end
            end
        end
        
        % Normalize 2nd level (capture contribution of each band to overall
        % power)
        % if norm
        if isfield(HRV.frequency,'ttlpwr')
            disp('Normalizing HRV power to overall power')
            if isfield(HRV.frequency,'ulf')
                HRV.frequency.ulf = mean(HRV.frequency.ulf) / HRV.frequency.ttlpwr;
            end
            if isfield(HRV.frequency,'vlf')
                HRV.frequency.vlf = mean(HRV.frequency.vlf) / HRV.frequency.ttlpwr;
            end
            if isfield(HRV.frequency,'lf')
                HRV.frequency.lf = mean(HRV.frequency.lf) / HRV.frequency.ttlpwr;
            end
            if isfield(HRV.frequency,'hf')
                HRV.frequency.hf = mean(HRV.frequency.hf) / HRV.frequency.ttlpwr;
            end
            if isfield(HRV.frequency,'lfhf')
                HRV.frequency.lfhf = round(mean(HRV.frequency.lf) / mean(HRV.frequency.hf) * 100)/100;
            end
        end
    % else
    %     warndlg('Sorry, HRV power could not be normalized to total power because ULF and VLF could not be estimated due to the short length of the cardiovascular time series.')
    %     warning('Sorry, HRV power could not be normalized to total power because ULF and VLF could not be estimated due to the short length of the cardiovascular time series.')
    end

    % remove empty cells
    PWR(cellfun(@isempty,PWR)) = []; 
    PWR_freqs(cellfun(@isempty,PWR_freqs)) = []; 

    % Merge spectra from each band and export
    HRV.frequency.pwr_freqs = [cat(1, PWR_freqs{:})];
    HRV.frequency.pwr = [cat(1, PWR{:})];
    HRV.frequency.bands = bands;
    
    % Average across time windows
    if isfield(HRV.frequency,'ulf') 
        HRV.frequency.ulf = round(mean(HRV.frequency.ulf,'omitnan'),2);
    end
    if isfield(HRV.frequency,'vlf') 
        HRV.frequency.vlf = round(mean(HRV.frequency.vlf,'omitnan'),2);
    end
    if isfield(HRV.frequency,'lf') 
        HRV.frequency.lf = round(mean(HRV.frequency.lf,'omitnan'),2);
    end
    if isfield(HRV.frequency,'hf') 
        HRV.frequency.hf = round(mean(HRV.frequency.hf,'omitnan'),2);
    end
    
    % if params.vis_outputs
    %     subplot(2,2,2+iWin); hold on;
    %     x = find(ulf_idx); y = pwr(ulf_idx);
    %     area(x,y,'FaceColor',[0.6350 0.0780 0.1840],'FaceAlpha',.7);
    %     x = find(vlf_idx); y = pwr(vlf_idx);
    %     area(x,y,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',.7)
    %     x = find(lf_idx); y = pwr(lf_idx);
    %     area(x,y,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',.7)
    %     x = find(hf_idx); y = pwr(hf_idx);
    %     area(x,y,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',.7)
    %     % xticks(1:8); xticklabels(reshape(fBounds,1,[]));
    %     xticks(1:30:length(f)); xticklabels(f(1:30:end));
    %     xlabel('Frequency (Hz)'); ylabel('Power (normalized)');
    %     legend('ULF', 'VLF', 'LF', 'HF')
    %     title(sprintf('Lomb-Scargle periodogram - Window %d', iWin))
    % end
end

%% Nonlinear domain
if params.hrv_nonlinear
    
    disp('Extracting HRV features in the nonlinear domain (Poincare, Sample entropy, Fractal dimension, PRSA)...')

    % Poincare
    SDSD = std(diff(NN));
    SDRR = std(NN);
    SD1 = (1 / sqrt(2)) * SDSD;     % measures the width of poincare cloud
    SD2 = sqrt((2 * SDRR^2) - (0.5 * SDSD^2));      % measures the length of the poincare cloud
    HRV.nonlinear.Poincare.SD1 = round(SD1*1000,3);      % in ms
    HRV.nonlinear.Poincare.SD2 = round(SD2*1000,3);      % in ms
    HRV.nonlinear.Poincare.SD1SD2 = round(SD1/SD2,3);

    % Sample entropy
    % m = 2;
    % r = .15;
    % HRV.nonlinear.SE = compute_se(NN,m,r);

    % Fuzzy entropy parameters
    m = 2;
    r = .15;
    tau = 1;
    n = 2;
    params.entropy_m = m;  % for exportation for users
    params.entropy_r = r;  % for exportation for users
    params.entropy_tau = tau;  % for exportation for users
    params.entropy_n = n;  % for exportation for users

    % Run Fuzzy entropy
    HRV.nonlinear.FE = compute_fe(NN, m, r, n, tau, false);

    % Multiscale fuzzy entropy (MFE)
    % coarseType = 'Standard deviation';
    % nScales = 30;
    % filtData = false;
    % [HRV.nonlinear.MFE, HRV.nonlinear.MFE_scales] = compute_mfe(NN, ...
    %     m, r, tau, coarseType, nScales, filtData, params.fs, n, params.gpu);
    % figure; area(HRV.MFE_scales, HRV.nonlinear.multiscale_fuzzy_entropy); axis tight

    % Fractal dimension
    HRV.nonlinear.FD = fractal_volatility(NN);
    
    % Phase rectified signal averaging (PRSA)
    fprintf('Computing phase rectified signal averaging (PRSA)... \n')
    thresh = 20;
    params.prsa_thresh = 20;  % for exportation for users
    lowAnchor = 1-thresh/100-0.0001; % lower limit for the AC anchor selection
    highAnchor = 1+thresh/100;      % The upper limit for the DC anchor selection
    drr_per = NN(2:end)./NN(1:end-1);
    ac_anchor = (drr_per > lowAnchor) & (drr_per <= .9999); % defines ac anchors, no changes greater than 5%
    dc_anchor = (drr_per > 1) & (drr_per <= highAnchor);
    ac_anchor(1) = false; % ignore 1st heartbeat
    dc_anchor(1) = false; % ignore 1st heartbeat
    HRV.nonlinear.PRSA_AC = round(mean(1000*NN(ac_anchor)),1);  % acceleration capacity (in ms)
    HRV.nonlinear.PRSA_DC = round(mean(1000*NN(dc_anchor)),1);  % deceleration capacity (in ms)
    
end


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

function NN_resamp = resample_NN(NN_times,NN,sf,interp_method)

% time index
ti = NN_times(1):1/sf:NN_times(end);     

% Resample with interpolation method of choice
switch interp_method
    case 'cub'
        NN_resamp = interp1(NN_times,NN,ti','spline')'; % cubic spline interpolation (default)
    case 'lin'
        NN_resamp = interp1(NN_times,NN,ti','linear')'; % linear interpolation
end
