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
    area(freqs(x),y,'FaceColor',"#A2142F",'FaceAlpha',.7);
else
    bandNames(idx) = [];
    bands(idx,:) = [];
end

idx = find(strcmp(bandNames,'VLF'));
if isfield(HRV.frequency, 'vlf')
    x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
    y = pwr(x);
    area(freqs(x),y,'FaceColor',"#D95319",'FaceAlpha',.7)
else
    bandNames(idx) = [];
    bands(idx,:) = [];
end

idx = find(strcmp(bandNames,'LF'));
if isfield(HRV.frequency, 'lf')
    x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
    y = pwr(x);
    area(freqs(x),y,'FaceColor',"#EDB120",'FaceAlpha',.7)
else
    bandNames(idx) = [];
    bands(idx,:) = [];
end

idx = find(strcmp(bandNames,'HF'));
if isfield(HRV.frequency, 'hf')
    x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
    y = pwr(x);
    area(freqs(x),y,'FaceColor',"#0072BD",'FaceAlpha',.7)
else
    bandNames(idx) = [];
    bands(idx,:) = [];
end
legend(bandNames)
% xticklabels(unique(reshape(bands,1,[]))); xtickangle(45); 
axis tight; box on
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
area(scales,mfe,'FaceColor',"#A2142F",'FaceAlpha',.7);
title('Multiscale fuzzy entropy - HRV'); xlabel('Scale factors'); ylabel('Entropy')
axis tight; box on; grid on

% PSD - EEG
subplot(2,2,2); hold on
bands = [0 3; 3 7; 7 13; 13 30; 30 max(freqs)];
pwr = trimmean(EEG.frequency.pwr,20,1); % 20% trimmed mean across channels
freqs = EEG.frequency.freqs(1,:);
% area(freqs,pwr,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',.7);

% delta
x = freqs >= bands(1,1) & freqs <= bands(1,2);
y = pwr(x);
area(freqs(x),y,'FaceColor',"#A2142F",'FaceAlpha',.7)

% theta
x = freqs >= bands(2,1) & freqs <= bands(2,2);
y = pwr(x);
area(freqs(x),y,'FaceColor',"#D95319",'FaceAlpha',.7)

% alpha
x = freqs >= bands(3,1) & freqs <= bands(3,2);
y = pwr(x);
area(freqs(x),y,'FaceColor',"#EDB120",'FaceAlpha',.7)

% beta
x = freqs >= bands(4,1) & freqs <= bands(4,2);
y = pwr(x);
area(freqs(x),y,'FaceColor',"#0072BD",'FaceAlpha',.7)

% gamma
x = freqs >= bands(5,1) & freqs <= bands(5,2);
y = pwr(x);
area(freqs(x),y,'FaceColor',"#4DBEEE",'FaceAlpha',.7)

% bandNames = {'Delta Theta Alpha Beta'};
% xticks = freqs;
% xticklabels(unique(reshape(bands,1,[])));
% xtickangle(45); axis tight; box on
title('Power spectral density - EEG'); legend({'delta' 'theta' 'alpha' 'beta' 'gamma'})
xlabel('Frequency (Hz)'); ylabel('Power (ms^2)'); axis tight;

% Multiscale fuzzy entropy (MFE) - EEG
subplot(2,2,4); %hold on
scales = EEG.nonlinear.MFE_scales(1,:); 
scaleBounds = EEG.nonlinear.MFE_scale_bounds(1,:); 
% mfe = mean(EEG.nonlinear.MFE,1); % mean across channels
mfe = trimmean(EEG.nonlinear.MFE,20,1); % 20% trimmed mean across channels
% mfe = EEG.nonlinear.MFE(31,:); % Pz
area(scales,mfe(end:-1:1),'FaceColor',"#A2142F",'FaceAlpha',.7);
axis tight; box on; grid on
if ~isempty(scaleBounds)
    xticks(scales); 
    xticklabels(scaleBounds(end:-1:1)); 
    xtickangle(45) %xticklabels({f}); %
else
    xlabel('Scale factors')
end
if max(mfe) < 1
    ylim([0 1])
end
title('Multiscale fuzzy entropy - EEG'); ylabel('Entropy')

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

%% Asymmetry 



%% coherence



% plot line between channels using normal topo and channel locations
% topoplot(results.diff, chanlocs, 'emarker2',{[3 17],'c','r'}); 

% % folder = 'C:\Users\IONSLAB\Desktop\channeling_matlab\plot_results';
% folder = 'C:\Users\Tracy\Documents\MATLAB\channeling_matlab';
% cd(fullfile(folder))
% 
% % Load head model and stat data
% hMdlPath = fileparts(which('head_modelColin27_5003_Standard-10-5-Cap339.mat'));
% atlas = load('-mat', fullfile(hMdlPath, 'head_modelColin27_5003_Standard-10-5-Cap339.mat'));
% areas = readtable('area_labels.csv');
% labels = {chanlocs.labels};
% 
% freqs = {'delta' 'theta' 'alpha' 'beta'};
% 
% % figure('color','w','position',[1353 904 1773 433])
% figure('color','w')
% for iFreq = 2:length(freqs)
%     disp(['Band: ' char(freqs(iFreq))])
% %     subplot(1,4,iFreq)
%     subplot(1,3,iFreq-1)
%     tmp = readtable(['pairwise_table_' freqs{iFreq} '.csv']);
% 
%     % Recompute pvals for theta with slightly more aggressive FDR-correction at 0.01 (cedric)
%     if strcmpi(freqs{iFreq}, 'theta')
%         disp([ 'Significant areas before new correction: ' num2str(sum(strcmp(tmp.Significant,'yes'))) ])
%         disp('New FDR-correction for theta at p=0.01: ')
%         idx = fdr_bh(tmp.pval,0.001,'pdep','yes'); 
%         tmp.Significant(idx==1) = {'yes'};
%         tmp.Significant(idx==0) = {'no'};
%     end
% 
%     array = zeros(68,68);
%     for iArea = 1:size(tmp,1)
%         if strcmpi(tmp.Significant{iArea}, 'yes')
%             array(tmp.area1(iArea), tmp.area2(iArea)) = tmp.diff(iArea);
%         end
%     end
%     iCol = 3; % 2 for long name; 3 for abreviations
% 
%     % Ordered by lobe
% %     plotconnectivity(array,'labels',areas(2:end,iCol)','brainimg','off', ...
% %         'threshold',0,'axis',gca,'labelsgroup',areas(2:end,4)'); 
% 
%     % Custom order to better approximate brain anatomy (front/back, left/right hemisphere)
%     plotconnectivity(array,'labels',areas(2:end,iCol)','brainimg','off', ...
%         'threshold',0,'axis',gca,'labelsgroup',areas(2:end,4)','reorder',areas(2:end,5));
% %     title(freqs(iFreq))
% end
% % legend(unique(areas(2:end,4)))
% 
% % print('-depsc', 'figure_connectivity.eps')
% % print('-djpeg', 'figure.jpg')
% % print(gcf,'figure.png');    %300 dpi
% % exportgraphics(gcf,'figure_300dpi.png','Resolution',300)  %300 dpi
% print(gcf,'figure.png','-dpng','-r300');                    %300 dpi
% print('figure_300dpi.tiff','-dtiff','-r300');               %300 dpi


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
