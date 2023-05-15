function plot_corrmatrix(C,labels)

load('corr_cmap.mat');

% format into a triangular matrix
C = tril(C,-1);                 % zero upper triangle
C(logical(eye(size(C)))) = 1;   % one in diagonal

% Set the  [min,max] of diameter where 1 consumes entire grid square
diamLim = [0.1, 1];

% Compute center of each circle
x = 1 : 1 : size(C,2); % x edges
y = 1 : 1 : size(C,1); % y edges
[xAll, yAll] = meshgrid(x,y);
xAll(C==0) = nan; % eliminate cordinates for zero correlations

% Scale colors from -1 to 1 (negative to positive correlations)
clrLim = [-1,1]; 
Cscaled = (C - clrLim(1))/range(clrLim); 
colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));

% Scale the size between 0-1
Cscaled = (abs(C) - 0)/1;
diamSize = Cscaled * range(diamLim) + diamLim(1);

% Create figure
fh = figure();
ax = axes(fh);
hold(ax,'on')
colormap(cmap) 
tickvalues = 1:size(C,2);
x = zeros(size(tickvalues));
text(x, tickvalues, labels,'HorizontalAlignment','right');
x(:) = size(C,1)+1;
text(tickvalues, x, labels,'HorizontalAlignment','right','Rotation',45);

% Create circles
theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
    diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:), ...
    'LineStyle','none'),1:numel(xAll));

% params
axis(ax,'equal'); axis(ax,'tight'); set(ax,'YDir','Reverse')
cb = colorbar(); ylabel(cb, 'Correlation coefficient','fontsize',12,'fontweight','bold')
clim(clrLim);
axis off;

set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
