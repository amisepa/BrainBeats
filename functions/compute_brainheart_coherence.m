% COMPUTE_BRAINHEART_COHERENCE - Brain-heart coherence and causal interactions.
%
% Fits an order-5 MVAR model to all channels (EEG + heart) by least squares
% (local fit_mvar), then computes the frequency-domain measures from its
% coefficients and residual covariance with fdMVAR_5order (Faes & Nollo 2011).
%
% Usage:
%   coherence = compute_brainheart_coherence(EEG, params)
%
% INPUTS:
%   EEG - EEGLAB structure containing EEG data (can be many channels) and
%       one cardiovascular channel (ECG or PPG). Both should be already
%       preprocessed to remove noise and artifacts.
%
%   params  - BrainBeats param structure. Only fields used by this function
%       are params.fs (sample rate), params.heart_channels (indicating the
%       name of the heart channel, e.g. {'ECG'}, for plots),
%       params.chanlocs (EEG channel locations for topo plots),
%       params.vis_outputs (to generate the plots or not).
%
% OUTPUTS:
%   coherence structure with fields freqs (1 x 4*fs, 0 to Nyquist in steps
%       of 1/8 Hz), channels (labels) and coh, pcoh, dc, pdc, gpdc (nChan x
%       nChan x nFreqs magnitudes, 0-1; element (i,j,:) is from channel j
%       to channel i):
%       - COH (coherence): linear coupling between two signals per frequency.
%       - PCOH (partial coherence): coupling after removing the influence of
%         all other channels.
%       - DC (directed coherence): total (direct + indirect) causal influence
%         of one signal on another.
%       - PDC (partial directed coherence): direct causal influence only.
%       - GPDC (generalized PDC): PDC weighted by the noise variances, so it
%         is not affected by differences in signal scale (e.g., ECG vs EEG).
%       DTF (directed transfer function) is also computed but not returned.
%
%   With params.vis_outputs, plots coh/pcoh/dc/pdc from EEG channels to the
%   heart channel (row = heart) up to 40 Hz, and scalp maps of their means in
%   delta (<= 3 Hz), theta (3-7), alpha (8-13) and beta (13-30 Hz).
%
% EXAMPLE USAGE:
%
%   params.fs = EEG.srate;
%   params.vis_outputs = 1;
%   params.heart_channels = {'ECG'};
%   params.chanlocs = EEG.chanlocs(~strcmpi({EEG.chanlocs.labels},'ECG'));  % EEG only
%   coherence = compute_brainheart_coherence(EEG,params)
%
% Copyright (C), Cedric Cannard, BrainBeats 2024
%
% PLEASE CITE THE FOLLOWING REFERENCE WHEN USING THIS CODE:
%   Faes & Nollo (2011). Multivariate Frequency Domain Analysis of Causal Interactions in Physiological Time Series. Biomedical Engineering, Trends in Electronics, Communications and Software.

function coherence = compute_brainheart_coherence(EEG,params)

disp('Computing brain-heart coherence measures...')
nfft = params.fs*4;  % number of frequency points from 0 to Nyquist (1/8 Hz resolution)

% Fit the order-5 MVAR model, then run the multivariate frequency-domain
% analysis of causal interactions on its coefficients (fdMVAR_5order takes
% model coefficients and noise covariance, not data)
[Am, Su] = fit_mvar(double(EEG.data), 5);
[dc,dtf,pdc,gpdc,~,coh,pcoh,~,~,~,~,f] = fdMVAR_5order(Am,Su,nfft,params.fs);

% Magnitudes (0-1) of the complex measures
if ~isreal(coh), coh = abs(coh); end
if ~isreal(pcoh),pcoh = abs(pcoh); end
if ~isreal(dc), dc = abs(dc); end
if ~isreal(pdc), pdc = abs(pdc); end
if ~isreal(gpdc), gpdc = abs(gpdc); end
if ~isreal(dtf), dtf = abs(dtf); end

% Outputs
coherence.freqs = f;
coherence.coh = coh;
coherence.pcoh = pcoh;
coherence.dc = dc;
coherence.pdc = pdc;
coherence.gpdc = gpdc;
coherence.channels = {EEG.chanlocs.labels};

% Plot
if params.vis_outputs

    disp("Plotting brain-heart coherence outputs...")
    figs0 = findall(groot, 'Type', 'figure');   % figures open before

    cardio_chan = strcmpi({EEG.chanlocs.labels},params.heart_channels);
    maxfreq = 40;

    % PLOT ALL FREQS AND CHANNELS FOR EACH MEASURE (heart row: EEG -> heart)
    figure('color','w')

    subplot(2,2,1) % COHERENCE
    fc = squeeze(coh(cardio_chan,~cardio_chan,f<=maxfreq));  % heart row, all frequencies up to maxfreq
    imagesc(f(f<maxfreq),1:size(fc,1),fc);  % EEG channels x frequencies
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Coherence','fontsize',12,'fontweight','bold','Rotation',270)
    Ylabels = {EEG.chanlocs(~cardio_chan).labels}; newticks = 1:2:length(Ylabels); newticks = unique(newticks);
    Ylabels  = Ylabels(newticks); set(gca,'YTick',newticks); set(gca,'YTickLabel', Ylabels,'FontWeight','normal');
    xlabel('Frequency (Hz)');

    subplot(2,2,2)  % PARTIAL COHERENCE
    fc = squeeze(pcoh(cardio_chan,~cardio_chan,f<=maxfreq));
    imagesc(f(f<maxfreq),1:size(fc,1),fc);  % EEG channels x frequencies
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Partial coherence','fontsize',12,'fontweight','bold','Rotation',270)
    Ylabels = {EEG.chanlocs(~cardio_chan).labels}; newticks = 1:2:length(Ylabels); newticks = unique(newticks);
    Ylabels  = Ylabels(newticks); set(gca,'YTick',newticks); set(gca,'YTickLabel', Ylabels,'FontWeight','normal');
    xlabel('Frequency (Hz)');

    subplot(2,2,3)  % DIRECTED COHERENCE
    fc = squeeze(dc(cardio_chan,~cardio_chan,f<=maxfreq));  % from EEG channels to heart
    imagesc(f(f<maxfreq),1:size(fc,1),fc);  % EEG channels x frequencies
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Directed coherence','fontsize',12,'fontweight','bold','Rotation',270)
    Ylabels = {EEG.chanlocs(~cardio_chan).labels}; newticks = 1:2:length(Ylabels); newticks = unique(newticks);
    Ylabels  = Ylabels(newticks); set(gca,'YTick',newticks); set(gca,'YTickLabel', Ylabels,'FontWeight','normal');
    xlabel('Frequency (Hz)');

    subplot(2,2,4)  % PARTIAL DIRECTED COHERENCE
    fc = squeeze(pdc(cardio_chan,~cardio_chan,f<=maxfreq));  % from EEG channels to heart
    imagesc(f(f<maxfreq),1:size(fc,1),fc);  % EEG channels x frequencies
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Partial directed coherence','fontsize',12,'fontweight','bold','Rotation',270)
    Ylabels = {EEG.chanlocs(~cardio_chan).labels}; newticks = 1:2:length(Ylabels); newticks = unique(newticks);
    Ylabels  = Ylabels(newticks); set(gca,'YTick',newticks); set(gca,'YTickLabel', Ylabels,'FontWeight','normal');
    xlabel('Frequency (Hz)');

    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

    % SCALP TOPOGRAPHY FOR EACH BAND: COHERENCE
    figure('color','w')
    subplot(2,2,1)  % delta
    fc = mean(coh(:,:,f<=3),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Delta','fontSize',12,'FontWeight','bold');
    subplot(2,2,2)  % theta
    fc = mean(coh(:,:,f>=3 & f<=7),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Coherence','fontsize',12,'fontweight','bold','Rotation',270)
    ylabel(cb,'Coherence','Rotation',270,'fontSize',12,'fontweight','bold')
    title('Theta','fontSize',12,'FontWeight','bold');
    subplot(2,2,3)  % alpha
    fc = mean(coh(:,:,f>=8 & f<=13),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Alpha','fontSize',12,'FontWeight','bold');
    subplot(2,2,4)  % beta
    fc = mean(coh(:,:,f>13 & f<=30),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Coherence','fontsize',12,'fontweight','bold','Rotation',270)
    ylabel(cb,'Coherence','Rotation',270,'fontSize',12,'fontweight','bold')
    title('Beta','fontSize',12,'FontWeight','bold');
    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

    % SCALP TOPOGRAPHY FOR EACH BAND: PARTIAL COHERENCE
    figure('color','w')
    subplot(2,2,1)  % delta
    fc = mean(pcoh(:,:,f<=3),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Partial coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Delta','fontSize',12,'FontWeight','bold');
    subplot(2,2,2)  % theta
    fc = mean(pcoh(:,:,f>=3 & f<=7),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Partial coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Theta','fontSize',12,'FontWeight','bold');
    subplot(2,2,3)  % alpha
    fc = mean(pcoh(:,:,f>=8 & f<=13),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Partial coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Alpha','fontSize',12,'FontWeight','bold');
    subplot(2,2,4)  % beta
    fc = mean(pcoh(:,:,f>13 & f<=30),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Partial coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Beta','fontSize',12,'FontWeight','bold');
    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

    % SCALP TOPOGRAPHY FOR EACH BAND: DIRECTED COHERENCE
    figure('color','w')
    subplot(2,2,1)  % delta
    fc = mean(dc(:,:,f<=3),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Directed coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Delta','fontSize',12,'FontWeight','bold');
    subplot(2,2,2)  % theta
    fc = mean(dc(:,:,f>=3 & f<=7),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Directed coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Theta','fontSize',12,'FontWeight','bold');
    subplot(2,2,3)  % alpha
    fc = mean(dc(:,:,f>=8 & f<=13),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Directed coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Alpha','fontSize',12,'FontWeight','bold');
    subplot(2,2,4)  % beta
    fc = mean(dc(:,:,f>13 & f<=30),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Directed coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Beta','fontSize',12,'FontWeight','bold');
    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

    % SCALP TOPOGRAPHY FOR EACH BAND: PARTIAL DIRECTED COHERENCE
    figure('color','w')
    subplot(2,2,1)  % delta
    fc = mean(pdc(:,:,f<=3),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Partial directed coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Delta','fontSize',12,'FontWeight','bold');
    subplot(2,2,2)  % theta
    fc = mean(pdc(:,:,f>=3 & f<=7),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Partial directed coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Theta','fontSize',12,'FontWeight','bold');
    subplot(2,2,3)  % alpha
    fc = mean(pdc(:,:,f>=8 & f<=13),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Partial directed coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Alpha','fontSize',12,'FontWeight','bold');
    subplot(2,2,4)  % beta
    fc = mean(pdc(:,:,f>13 & f<=30),3,'omitnan');  % mean for the band
    plot_topo(fc(cardio_chan,~cardio_chan), params.chanlocs, 1, 'psd');  % cardio in row, EEG in columns
    set(gca,'CLim',[0 1]); cb = colorbar; ylabel(cb,'Partial directed coherence','fontsize',12,'fontweight','bold','Rotation',270)
    title('Beta','fontSize',12,'FontWeight','bold');
    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

    figs = findall(groot, 'Type', 'figure');
    finish_figure(figs(~ismember(figs, figs0)))
end

disp("Done computing Coherence measures.")


%% Subfunction
function [Am, Su] = fit_mvar(X, p)
% Least-squares identification of a strictly causal MVAR model of order p
% (X: channels x samples), as idMVAR in Faes & Nollo (2011).
%   Am = [A(1) ... A(p)] (M x pM coefficients), Su = M x M residual covariance
[M, N] = size(X);
X = X - mean(X,2);
Y = X(:, p+1:N);
Z = zeros(p*M, N-p);
for k = 1:p
    Z((k-1)*M+1:k*M, :) = X(:, p+1-k:N-k);
end
Am = (Y*Z') / (Z*Z');
U  = Y - Am*Z;
Su = (U*U') / (N-p);
