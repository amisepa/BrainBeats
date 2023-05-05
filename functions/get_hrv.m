function HRV = get_hrv(NN, NN_times, params)


% Time domain
if params.hrv_time
    % TimeMetrics = EvalTimeDomainHRVstats(NN,tNN,sqi,HRVparams,tWin);
    HRV.time = compute_time(NN, params);
end

% Frequency domain
if params.hrv_frequency
    % ADD: CHECK IF BANDS CAN BE COMPUTED WITH FILE LENGTH AVAILABLE
    HRV.frequency = compute_frequency(NN, NN_times, params);
end

% Poincare
if params.hrv_nonlinear
    SDSD = std(diff(NN));
    SDRR = std(NN);
    SD1 = (1 / sqrt(2)) * SDSD;     % measures the width of poincare cloud
    SD2 = sqrt((2 * SDRR^2) - (0.5 * SDSD^2));      % measures the length of the poincare cloud
    HRV.Poincare.SD1 = SD1*1000;      % in ms
    HRV.Poincare.SD2 = SD2*1000;      % in ms
    HRV.Poincare.SD1SD2 = SD1/SD2;

    % % Entropy FIXME: ADD CHECK THAT get_entropy plugin is installed in EEGLAB
    % % Parameters
    % tau = 1;
    % m = 2;
    % coarseType = 'Standard deviation';
    % nScales = 30;
    % r = .15;
    % n = 2;
    % filtData = false;
    % useGPU = false;
    % 
    % HRV.appEn = compute_ae(zscore(NN), m, r);
    % HRV.sampEn = compute_se(zscore(NN), m, r, tau );
    % HRV.fuzzEn = compute_fe(zscore(NN), m, r, n, tau);
    % [HRV.MSE, HRV.MSE_scales] = compute_mse(zscore(NN), ...
    %     m, r, tau, coarseType, nScales, filtData, srate); % ADD remove NaNs
    % [HRV.MFE, HRV.MSE_scales] = compute_mfe(zscore(NN), ...
    %     m, r, tau, coarseType, nScales, filtData, srate, n, useGPU); % ADD remove NaNs


% Phase rectified signal averaging (PRSA) (FIXME)
% if params.hrv_prsa
%
%     % parameters
%     thresh = 20;
%     min_anch = 20;
%     scale = 2;
%
%     lowAnchor = 1-thresh/100-.0001; % lower limit for the AC anchor selection
%     highAnchor = 1+thresh/100;      % The upper limit for the DC anchor selection
%     drr_per = NN(2:end)./NN(1:end-1);
%     ac_anchor = (drr_per > lowAnchor) & (drr_per <= .9999); % defines ac anchors, no changes greater than 5%
%     dc_anchor = (drr_per > 1) & (drr_per <= highAnchor);
%     ac_anchor = [0; ac_anchor(:)];
%     dc_anchor = [0; dc_anchor(:)];
%     ac_anchor(1:winSize) = 0;                                        % sets first window to 0
%     ac_anchor(length(ac_anchor)-winSize+1:length(ac_anchor)) = 0;    % sets last window to 0
%     dc_anchor(1:winSize) = 0;                                        % sets first window to 0
%     dc_anchor(length(dc_anchor)-winSize+1:length(dc_anchor)) = 0;    % sets last window to 0
%     ac_ind = find(ac_anchor);
%     dc_ind = find(dc_anchor);
%     for i = 1:length(dc_ind)
%         dcm(i,:) = 1000*NN(dc_ind(i)-winSize:dc_ind(i)+winSize-1); % in ms
%     end
%     for i = 1:length(ac_ind)
%         acm(i,:) = 1000*NN(ac_ind(i)-winSize:ac_ind(i)+winSize-1); % convert from sec to msec
%     end
%     prsa_ac = mean(acm,1);
%     prsa_dc = mean(dcm,1);
%
%     % For unusual NN intervals
%     if ~isnan(sum(prsa_dc)) && numel(dc_ind) >= thresh
%         dc = (sum(prsa_dc(winSize+1:winSize+scale)) - sum(prsa_dc(winSize-(scale-1):winSize))) ./ (2*scale);
%         dc_results = dc; % assign output of window
%     end
%     if ~isnan(sum(prsa_ac)) && numel(ac_ind) >= min_anch
%         ac = (sum(prsa_ac(winSize+1:winSize+scale)) - sum(prsa_ac(winSize-(scale-1):winSize))) ./ (2*scale);
%         ac_results = ac; % assign output of window
%     end
%
%     HRV.prsa_ac = ac_results; % acceleration capacity
%     HRV.prsa_dc = dc_results; % deceleration capacity
end

% Visualization of HRV outputs
if params.vis
    figure('color','w');
    
    % Power spectra
    % subplot(2,1,1); 
    hold on
    % if nWin == 1
        pwr = HRV.frequency.pwr;
    % else
        % pwr = mean([HRV.frequency.pwr],1);
    % end
    x = find(HRV.frequency(1).ulf_idx);
    y = pwr(HRV.frequency(1).ulf_idx);
    area(x,y,'FaceColor',[0.6350 0.0780 0.1840],'FaceAlpha',.7);
    my_ticks = [x(1) x(end)];
    x = find(HRV.frequency(1).vlf_idx);
    y = pwr(HRV.frequency(1).vlf_idx);
    area(x,y,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',.7)
    my_ticks(3) = x(end);
    x = find(HRV.frequency(1).lf_idx);
    y = pwr(HRV.frequency(1).lf_idx);
    area(x,y,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',.7)
    my_ticks(4) = x(end);
    x = find(HRV.frequency(1).hf_idx);
    y = pwr(HRV.frequency(1).hf_idx);
    area(x,y,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',.7)
    my_ticks(5) = x(end);
    xticks(my_ticks);
    xticklabels(unique(reshape(HRV.frequency(1).bands,1,[])));
    xtickangle(45); axis tight; box on
    xlabel('Frequency (Hz)');
    if params.hrv_norm
        ylabel('Power (ms^2 normalized)');
    else
        ylabel('Power (ms^2)');
    end
    legend('ULF', 'VLF', 'LF', 'HF')
    title(sprintf('Lomb-Scargle periodogram'))

    % Multiscale fuzzy entropy (MFE)
    % subplot(2,1,1); hold on

    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
    % set(gca,'FontSize',12,'layer','top','fontweight','bold');

end
