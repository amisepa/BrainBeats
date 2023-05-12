function plot_features(Features)

HRV = Features.HRV;
eeg_features = Features.EEG;

figure('color','w');

% PSD - HRV
subplot(2,2,1); hold on
pwr = HRV.frequency.pwr; % power already averaged across windows
freqs = HRV.frequency.pwr_freqs;
bandNames = {'ULF' 'VLF' 'LF' 'HF'};
bands = HRV.frequency.bands;

idx = find(strcmp(bandNames,'ULF'));
if isfield(HRV.frequency, 'ulf')
    x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
    y = pwr(x);
    area(find(x),y,'FaceColor',[0.6350 0.0780 0.1840],'FaceAlpha',.7);
else
    bandNames(idx) = [];
    bands(idx,:) = [];
end

idx = find(strcmp(bandNames,'VLF'));
if isfield(HRV.frequency, 'vlf')
    x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
    y = pwr(x);
    area(find(x),y,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',.7)
else
    bandNames(idx) = [];
    bands(idx,:) = [];
end

idx = find(strcmp(bandNames,'LF'));
if isfield(HRV.frequency, 'lf')
    x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
    y = pwr(x);
    area(find(x),y,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',.7)
else
    bandNames(idx) = [];
    bands(idx,:) = [];
end

idx = find(strcmp(bandNames,'HF'));
if isfield(HRV.frequency, 'hf')
    x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
    y = pwr(x);
    area(find(x),y,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',.7)
else
    bandNames(idx) = [];
    bands(idx,:) = [];
end
legend(bandNames)
xticklabels(unique(reshape(bands,1,[])));
xtickangle(45); axis tight; box on
xlabel('Frequency (Hz)');
if params.hrv_norm
    ylabel('Power (ms^2 normalized)');
else
    ylabel('Power (ms^2)');
end
title(sprintf('Power spectral density - HRV'))

% Multiscale fuzzy entropy (MFE) - HRV
subplot(2,2,3); hold on
mfe = HRV.nonlinear.MFE;
scales = HRV.nonlinear.MFE_scales;
area(scales,mfe,'FaceColor',[0.6350 0.0780 0.1840],'FaceAlpha',.7);
title('Multiscale fuzzy entropy - HRV'); xlabel('Time scales'); ylabel('Entropy')
axis tight;

% PSD - EEG
% load('C:\Users\Tracy\Desktop\eeg_features.mat')
subplot(2,2,2); hold on
pwr = trimmean(eeg_features.frequency.pwr,20,1); % mean across channels
freqs = eeg_features.frequency.freqs(1,:);
area(freqs,pwr,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',.7);
title('Power spectral density - EEG'); xlabel('Frequency (Hz)');
ylabel('Power (ms^2)'); axis tight;

% Multiscale fuzzy entropy (MFE) - EEG
subplot(2,2,4); hold on
mfe = trimmean(eeg_features.nonlinear.MFE,20,1);
scales = eeg_features.nonlinear.MFE_scales(1,:);
area(scales,mfe,'FaceColor',[0.6350 0.0780 0.1840],'FaceAlpha',.7);
title('Multiscale fuzzy entropy - HRV'); xlabel('Time scales'); ylabel('Entropy')
axis tight;

set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
% set(gca,'FontSize',12,'layer','top','fontweight','bold');
