%% Compute HRV measures on NN time series. 
% 1) Time domain: NN stats, SDNN, RMSSD, pNN50
% 
% 2) Frequency domain: ULF, VLF, LF, HF, LF/HF ratio, TTLPWR.
%   Options: PSD can be calculated using the Lomb-Scargle periodogram (default), 
%   pwelch, FFT, or Burg methods. Lomb-Scargle periodogram does not require 
%   interpolation or resampling of the data (contrary to welch or FFT), thus 
%   preserving the original information. 
%   Normalization can be turned ON to better deal with non-uniformly sampled
%   or irregularly sampled data, which is common in HRV analysis. It is
%   applied during Lomb-scargle periodogram estimation by scaling the total
%   power with variance in the time series, to make results more comparable 
%   across different datasets or subjects. It is also applied in a 2nd
%   step by dividing each band power by the total power, to provide a more 
%   intuitive measure of the relative contribution of each frequency 
%   component to the overall power.
%
% 3) Nonlinear domain: Poincare, fuzzy entropy, multiscale fuzzy entropy,
%   and PRSA.
%
% Following recommendations by the Task Force of the European Society of 
% Cardiology and the North American Society of Pacing and Electrophysiology
% (1996), minimum data length for each band is 5-10 cycles.
%   - ULF: at least 24 hours
%   - VLF: at least 28 minutes
%   - LF: at least 125 s
%   - HF: at least 34 s
% Different window lengths are implemented to maximize trade-off between 
% time resolution and frequency resolution for each band.
% 
% Currently, HRV metrics are computed over the whole time-series for faster 
% computation and allow computation of ULF and VLF as much as possible
% (requiring long recordings to be accurate). Sliding time windows will be 
% implemented in the future for assessing changes over time, when enough 
% data are provided. 
% 
% Cedric Cannard, 2023

function HRV = get_hrv_features(NN, NN_times, params)

% FIXME: add sliding windows for users wanting to look at how it changes over time.

%% Time domain
if params.hrv_time

    disp('Extracting HRV features in the time domain...')

    % if NN_times(end) < 300
    %     warning('File length is too short to reliably estimate HRV time-domain metrics. At least 5 minutes of data are recommended (see Shaffer and Ginsberg, 2017).')
    % end

    % NN statistics
    HRV.time.NN_mean = mean(NN.*1000);      % in ms
    HRV.time.NN_median = median(NN.*1000);  % in ms
    HRV.time.NN_mode = mode(NN.*1000);      % in ms
    HRV.time.NN_var = var(NN.*1000);        % in ms
    HRV.time.NN_skew = skewness(NN);
    HRV.time.NN_kurt = kurtosis(NN);
    HRV.time.NN_iqr = iqr(NN.*1000);        % in ms

    % SDNN
    HRV.time.SDNN = std(NN.*1000);   % in ms

    % RMSSD (sqrt of the mean squared time diff between heartbeats)
    HRV.time.RMSSD = sqrt(mean(diff(NN.*1000).^2));  % in ms

    % pNN50 (fraction of differences larger than alpha = 50)
    alpha = 50;
    HRV.time.pNN50 = sum( abs(diff(NN)) >= alpha/1000 )/length(diff(NN));

end

%% Frequency domain
if params.hrv_frequency

    disp('Extracting HRV features in the frequency domain...')

    % HRV frequency bands (ULF; VLF; LF; HF)
    bands = [ 0 .003; 0.003 .04; .04 .15; 0.15 0.40 ];
    bandNames = {'ULF' 'VLF' 'LF' 'HF'};

    % Minimum data length requirements for each band
    % minULF = 86400;    % 24 hours
    % minVLF = 5/0.003;    % 5 cycles/0.03 hz  (in s)
    % minLF = 5/0.04;     % 5 cycles/0.04 hz  (in s)
    % minHF = 5/0.15;     % 5 cycles/0.15 hz  (in s)
    minLength = ceil([ 86400  5/0.003 5/0.04 5/0.15 ]);
    
    for iBand = 1:size(bands,1)

        if NN_times(end) >= minLength(iBand)
    
            % Determine best sliding window length and indices
            winLength = minLength(iBand);
            stepSize = winLength * (1 - params.hrv_overlap);
            nWindows = floor((NN_times(end) - winLength) / stepSize) + 1;
    
            % Compute PSD on each sliding window
            for iWin = 1:nWindows
                start_idx = (iWin - 1) * stepSize + 1;
                end_idx = start_idx + winLength - 1;
                win_idx = NN_times >= start_idx & NN_times <= end_idx;

                % Frequency vector for this band
                nfft = 1024;    % use this instead? 2^nextpow2(length(NN(win_idx)))
                fvec = bands(iBand,1):1/nfft:bands(iBand,2);
                                
                % Lomb-Scargle Periodogram (no resampling required; best method)
                if strcmp(params.hrv_spec, 'Lomb-Scargle periodogram')
                    if params.hrv_norm
                        [pwr,freqs] = plomb(NN(win_idx),NN_times(win_idx),fvec,'normalized'); 
                        fprintf('Computing normalized Lomb-Scargle periodogram on the NN series. \n')
                    else
                        [pwr,freqs] = plomb(NN(win_idx),NN_times(win_idx),fvec,'psd'); 
                        fprintf('Computing standard Lomb-Scargle periodogram on the NN series. \n')
                    end 
                
                % pwelch or FFT (require resampling)
                else
                    resamp_freq = 7;

                    % Resample
                    NN_resamp = resample_NN(NN_times(win_idx),NN(win_idx),resamp_freq,'cub'); % resample RR
            
                    % Pwelch
                    if strcmp(params.hrv_spec, 'pwelch')
                        [pwr,freqs] = pwelch(NN_resamp,[],[],2^nextpow2(length(NN_resamp)),resamp_freq);
                        
                    % FFT
                    elseif strcmp(params.hrv_spec, 'fft')
                        pwr = fft(NN_resamp).*conj(fft(NN_resamp))/length(NN_resamp);
                        freqs = resamp_freq*(0:floor(length(NN_resamp)/2)+1)/length(NN_resamp);
                        pwr = pwr(1:length(freqs));
            
                    else
                        error("params.hrv_spec can only be 'Lomb-Scargle periodogram', 'pwelch', or 'fft'. ")
                    end
                end
                
                % Freq index
                freq_idx = bands(iBand,1) <= freqs & freqs <= bands(iBand,2);
                freq_res = freqs(2)-freqs(1); % resolution

                % Power for each band in ms^2
                if iBand == 1
                    % HRV.frequency.ulf_idx = freq_idx;
                    HRV.frequency.ulf(iWin,:) = sum(pwr(freq_idx)*freq_res) * 1e6;    % ULF
                elseif iBand == 2
                    % HRV.frequency.vlf_idx = freq_idx;
                    HRV.frequency.vlf(iWin,:) = sum(pwr(freq_idx)*freq_res) * 1e6;    % VLF
                elseif iBand == 3
                    % HRV.frequency.lf_idx = freq_idx;
                    HRV.frequency.lf(iWin,:) = sum(pwr(freq_idx)*freq_res) * 1e6;    % LF
                elseif iBand == 4
                    % HRV.frequency.hf_idx = freq_idx;
                    HRV.frequency.hf(iWin,:) = sum(pwr(freq_idx)*freq_res) * 1e6;    % HF
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
    if isfield(HRV.frequency,'ulf') && isfield(HRV.frequency,'vlf') ...
            && isfield(HRV.frequency,'lf') && isfield(HRV.frequency,'hf')

        HRV.frequency.ttlpwr = sum([mean(HRV.frequency.ulf) mean(HRV.frequency.vlf) ...
            mean(HRV.frequency.lf) mean(HRV.frequency.hf)]);
        
        % Normalize 2nd level (capture contribution of each band to overall
        % power)
        if params.norm
            disp('Normalizing power')
            HRV.frequency.ulf = mean(HRV.frequency.ulf) / HRV.frequency.ttlpwr;
            HRV.frequency.vlf = mean(HRV.frequency.vlf) / HRV.frequency.ttlpwr;
            HRV.frequency.lf = mean(HRV.frequency.lf) / HRV.frequency.ttlpwr;
            HRV.frequency.hf = mean(HRV.frequency.hf) / HRV.frequency.ttlpwr;
            HRV.frequency.lfhf = round(mean(HRV.frequency.lf) / mean(HRV.frequency.hf) * 100)/100;
        end
    end

    % Merge spectra from each band and export
    PWR(cellfun(@isempty,PWR)) = []; PWR_freqs(cellfun(@isempty,PWR_freqs)) = []; 
    HRV.frequency.pwr_freqs = [cat(1, PWR_freqs{:})];
    HRV.frequency.pwr = [cat(1, PWR{:})];
    HRV.frequency.bands = bands;

    % if vis
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

% Nonlinear domain
if params.hrv_nonlinear

    % Poincare
    SDSD = std(diff(NN));
    SDRR = std(NN);
    SD1 = (1 / sqrt(2)) * SDSD;     % measures the width of poincare cloud
    SD2 = sqrt((2 * SDRR^2) - (0.5 * SDSD^2));      % measures the length of the poincare cloud
    HRV.nonlinear.Poincare.SD1 = SD1*1000;      % in ms
    HRV.nonlinear.Poincare.SD2 = SD2*1000;      % in ms
    HRV.nonlinear.Poincare.SD1SD2 = SD1/SD2;

    % Entropy (FIXME: install get_entropy plugin if not already installed)
    tau = 1;
    m = 2;
    coarseType = 'Standard deviation';
    nScales = 20;
    r = .15;
    n = 2;
    filtData = false;
    useGPU = false;
    % Fuzzy entropy
    HRV.nonlinear.FE = compute_fe(NN, m, r, n, tau,useGPU);
    % Multiscale fuzzy entropy
    [HRV.nonlinear.MFE, HRV.nonlinear.MFE_scales] = compute_mfe(NN, ...
        m, r, tau, coarseType, nScales, filtData, params.fs, n, useGPU);
    % figure; area(HRV.MFE_scales, HRV.nonlinear.multiscale_fuzzy_entropy); axis tight
    


    % Phase rectified signal averaging (PRSA) (FIXME)
    thresh = 20;
    lowAnchor = 1-thresh/100-.0001; % lower limit for the AC anchor selection
    highAnchor = 1+thresh/100;      % The upper limit for the DC anchor selection
    drr_per = NN(2:end)./NN(1:end-1);
    ac_anchor = (drr_per > lowAnchor) & (drr_per <= .9999); % defines ac anchors, no changes greater than 5%
    dc_anchor = (drr_per > 1) & (drr_per <= highAnchor);
    ac_anchor(1) = false; % ignore 1st heartbeat
    dc_anchor(1) = false; % ignore 1st heartbeat
    HRV.nonlinear.PRSA_AC = mean(1000*NN(ac_anchor));  % acceleration capacity (in ms)
    HRV.nonlinear.PRSA_DC = mean(1000*NN(dc_anchor));  % deceleration capacity (in ms)

end

end

%% Resample NN intervals to compute PSD with pwelch or FFT
%
%   INPUT:
%       win_idx         - start and end time of the NN interval
%       NN              - vector of NN intervals to be resampled
%       sf              - 7
%       interp_method   - 'cub' 
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
end