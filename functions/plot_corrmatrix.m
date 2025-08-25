function plot_corrmatrix(C, labels)
% plot_corrmatrix  Visualize a correlation matrix with colored circles
%
% Syntax
%   plot_corrmatrix(C, labels)
%
% Description
%   Draws a lower triangular correlation plot where circle color encodes
%   the sign and magnitude of correlation and circle size encodes the
%   absolute value. The upper triangle is hidden to reduce clutter.
%
% Inputs
%   C       Square correlation matrix of size [N x N], values in [-1, 1]
%   labels  Cell array or string array of length N with variable names
%
% Notes
%   1) The function loads a colormap from corr_cmap.mat and expects a
%      variable named cmap of size [K x 3] with RGB values in [0, 1].
%   2) Only the lower triangle is rendered. Diagonal entries are not drawn.
%   3) If C and labels sizes do not match, the function throws an error.
%
% Example: 
%   Nvars = 8;
%   Nsamp = 100;
%   X = randn(Nsamp, Nvars);
%   X(:,3) = 0.6*X(:,1) + 0.2*X(:,2) + 0.6*randn(Nsamp,1); % add weak structure
%   X(:,5) = -0.5*X(:,2) + 0.2*X(:,4) + 0.6*randn(Nsamp,1); % add weak structure
%   C = corrcoef(X);
%   labels = arrayfun(@(k)sprintf('Var%d',k), 1:Nvars, 'UniformOutput', false);
%   plot_corrmatrix(C, labels)
% 
% Cedric Cannard, 2023

if size(C,1) ~= length(labels)
    error("Labels and correlation matrix must have the same number of variables: %g", size(C,1))
end

load('corr_cmap.mat');  % loads variable 'cmap'

% format into a triangular matrix
% C = tril(C,-1);  % zero upper triangle (+ remove diagonal)
C = tril(C,0);    % zero upper triangle, preserving diagonal

% Set the [min,max] of diameter where 1 consumes entire grid square
diamLim = [0.1, 1];

% Compute center of each circle
x = 1 : 1 : size(C,2); % x edges
y = 1 : 1 : size(C,1); % y edges
[xAll, yAll] = meshgrid(x,y);
xAll(C==0) = nan; % eliminate coordinates for zero correlations

% Scale colors from -1 to 1 (negative to positive correlations)
clrLim = [-1,1];
Cscaled = (C - clrLim(1))/range(clrLim);
colIdx = discretize(Cscaled, linspace(0,1,size(cmap,1)));

% Scale the size between 0 and 1
Cscaled = abs(C);
diamSize = Cscaled * range(diamLim) + diamLim(1);

% Create figure
fh = figure();
ax = axes(fh);
hold(ax,'on')
colormap(cmap)
tickvalues = 1:size(C,2);
x = zeros(size(tickvalues));
text(x, tickvalues, labels,'HorizontalAlignment','right','fontSize',12,'fontweight','normal');
x(:) = size(C,1)+1;
text(tickvalues, x, labels,'HorizontalAlignment','right','Rotation',90,'fontSize',12,'fontweight','normal');

% Create circles
theta = linspace(0,2*pi,50); % the smaller, the less memory
arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
    diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:), ...
    'LineStyle','none'), 1:numel(xAll));

% Figure params
set(ax,'YDir','Reverse')
cb = colorbar();
ylabel(cb, 'Correlation coefficient','fontsize',12,'fontweight','bold')
clim(clrLim);
axis off;
set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold');
set(gcf,'Name','Correlation matrix','color','w','Toolbar','none','Menu','none','NumberTitle','Off')

end
