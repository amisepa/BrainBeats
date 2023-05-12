function plot_features(Features)

HRV = Features.HRV;
eeg_features = Features.EEG;


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

%% CORRELATION PLOT 1 (https://www.mathworks.com/help/econ/corrplot.html)

% pairwise Pearson's correlations and corresponding p-values for testing 
% the null hypothesis of no correlation against the right-tailed alternative 
% that the correlations are greater than zero. 

% load Data_Canada
% [R,PValue] = corrplot(DataTable,Tail="right");
% PValue

%% CORRELATION PLOT 2 (https://www.mathworks.com/matlabcentral/answers/699755-fancy-correlation-plots-in-matlab)

% % Produce the input lower triangular matrix data
% C = -1 + 2.*rand(12,12);
% C = tril(C,-1);
% C(logical(eye(size(C)))) = 1;
% 
% % Set [min,max] value of C to scale colors
% clrLim = [-1,1];
% 
% % load('CorrColormap.mat') % Uncomment for custom CorrColormap
% 
% % Set the  [min,max] of diameter where 1 consumes entire grid square
% diamLim = [0.1, 1];
% myLabel = {'ICA','Elev','Pr','Rmax','Rmin','Srad','Wspd','Tmin','Tmax','VPD','ET_o','AW'};
% 
% % Compute center of each circle (this assumes the x and y values were not 
% % centered in imagesc()
% x = 1 : 1 : size(C,2); % x edges
% y = 1 : 1 : size(C,1); % y edges
% [xAll, yAll] = meshgrid(x,y);
% xAll(C==0)=nan; % eliminate cordinates for zero correlations
% % Set color of each rectangle
% % Set color scale
% cmap = jet(256);
% % cmap = CorrColormap; % Uncomment for CorrColormap
% Cscaled = (C - clrLim(1))/range(clrLim); % always [0:1]
% colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% % Set size of each circle
% % Scale the size between [0 1]
% Cscaled = (abs(C) - 0)/1;
% diamSize = Cscaled * range(diamLim) + diamLim(1);
% % Create figure
% fh = figure();
% ax = axes(fh);
% hold(ax,'on')
% colormap(ax,'jet');
% % colormap(CorrColormap) %Uncomment for CorrColormap
% tickvalues = 1:length(C);
% x = zeros(size(tickvalues));
% text(x, tickvalues, myLabel, 'HorizontalAlignment', 'right');
% x(:) = length(C)+1;
% text(tickvalues, x, myLabel, 'HorizontalAlignment', 'right','Rotation',90);
% % Create circles
% theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
% h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
%     diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:),'LineStyle','none'),1:numel(xAll));
% axis(ax,'equal')
% axis(ax,'tight')
% set(ax,'YDir','Reverse')
% colorbar()
% clim(clrLim);
% axis off
