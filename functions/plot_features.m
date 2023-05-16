function plot_features(Features,params)

HRV = Features.HRV;
EEG = Features.EEG;


%% PSD AND MFE PLOT

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
subplot(2,2,2); hold on
pwr = trimmean(EEG.frequency.pwr,20,1); % 20% trimmed mean across channels
freqs = EEG.frequency.freqs(1,:);
area(freqs,pwr,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',.7);
title('Power spectral density - EEG'); xlabel('Frequency (Hz)');
ylabel('Power (ms^2)'); axis tight;

% Multiscale fuzzy entropy (MFE) - EEG
subplot(2,2,4); hold on
scales = EEG.nonlinear.MFE_scales(1,:); 
mfe = trimmean(EEG.nonlinear.MFE,20,1); % 20% trimmed mean across channels
% mfe = EEG.nonlinear.MFE(31,:); % elec with more variation
area(scales,mfe,'FaceColor',[0.6350 0.0780 0.1840],'FaceAlpha',.7);
title('Multiscale fuzzy entropy - HRV'); xlabel('Time scales'); ylabel('Entropy')
axis tight;

set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
% set(gca,'FontSize',11,'layer','top','fontweight','bold');

%% Scalp topos 2D

mode = 1;
figure('color','w')
subplot(3,2,1)
plot_topo(gather(mean(EEG.frequency.delta,2)),params.chanlocs,mode,'psd');
title('Delta power'); %colorbar off; 
subplot(3,2,2)
plot_topo(gather(mean(EEG.frequency.theta,2)),params.chanlocs,mode,'psd');
title('Theta power'); 
subplot(3,2,3)
plot_topo(gather(mean(EEG.frequency.alpha,2)),params.chanlocs,mode,'psd');
title('Alpha power'); 
subplot(3,2,4)
plot_topo(gather(mean(EEG.frequency.beta,2)),params.chanlocs,mode,'psd');
title('Beta power'); 
subplot(3,2,5)
plot_topo(gather(EEG.frequency.IAF),params.chanlocs,mode,'IAF');
title('Individual alpha frequency (IAF)'); %colorbar off; 
subplot(3,2,6)
plot_topo(gather(EEG.nonlinear.FE),params.chanlocs,mode,'entropy');
title('Fuzzy entropy'); %colorbar off; 

set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');


%% Scalp topos 3D - PSD and MFE

% mode = 2;
% dataToPlot{1,1} = gather(EEG.frequency.pwr_dB);
% dataToPlot{2,1} = gather(EEG.frequency.freqs);
% dataToPlot{1,2} = gather(EEG.nonlinear.MFE);
% dataToPlot{2,2} = gather(EEG.nonlinear.MFE_scales);
% figure('color','w')
% plot_topo(dataToPlot,params.chanlocs,mode,'psd');
% title('PSD and MFE'); %colorbar off; 
% set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');

%% Asymmetry and coherence





%% CORRELATION PLOT 

% labels = {'ICA','Elev','Pr','Rmax','Rmin','Srad','Wspd','Tmin','Tmax','VPD','ET_o','AW'};
% C = -1 + 2.*rand(12,12);      % produce fake data

% load hospital
% X = [hospital.Weight hospital.BloodPressure];
% R = corrcoef(X);


% X(1,2) = [HRV.time.NN_mean];
% X(2,1) = [HRV.time.NN_mean];

% plot_corrmatrix(C,labels)





%% CORRELATION PLOT 2 (requires econometrics toolbox)
% (https://www.mathworks.com/help/econ/corrplot.html)

% pairwise Pearson's correlations and corresponding p-values for testing 
% the null hypothesis of no correlation against the right-tailed alternative 
% that the correlations are greater than zero. 

% load Data_Canada
% [R,PValue] = corrplot(DataTable,Tail="right");
% PValue
