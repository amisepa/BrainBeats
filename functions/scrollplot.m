function scrollplot(varargin)

% Extract inputs
toPlot1 = varargin{1};      % Data to plot
scrollOrient = varargin{2}; % 'X' for horizontal and 'Y' for vertical
toDisplay = varargin{3};    % Window size (in s)

% Additional data to plot on top (e.g., markers)
toPlot2 = [];
toPlot3 = [];

if nargin >= 4
    toPlot2 = varargin{4};
end
if nargin == 5
    toPlot3 = varargin{5};
end

% Convert window size to percentage
toDisplay = toDisplay / max(toPlot1{1});

% Plot primary data
p = plot(toPlot1{:}, 'linewidth', 1);
ax = p.Parent;

% Set initial axes limits to a portion of the total plotted X data
XX = p.XData;
if ~isempty(XX)
    YY = p.YData;
    if ~all(isnan(XX))
        totxlim = [min(XX, [], 'omitnan'), max(XX, [], 'omitnan')];
        xlim = [totxlim(1), totxlim(1) + toDisplay * diff(totxlim)]; % Arbitrarily show first segment
    end
    if ~all(isnan(YY))
        % ylim = [-3*std(YY, [], 'omitnan'), 3*std(YY, [], 'omitnan')];
        % ylim = [min(YY)*1.1, max(YY)*1.1];
        ylim = [-std(YY,[], 'omitnan'), std(YY, [],'omitnan')];
    end
end
set(ax, 'xlim', xlim, 'ylim', ylim);

% Plot additional data (e.g., markers)
hold on
if ~isempty(toPlot2)
    plot(toPlot2{:});
end
if ~isempty(toPlot3)
    plot(toPlot3{:});
end

% Initialize scrolling behavior
superscroll_obj = superscroll(ax, scrollOrient);

% Build scrollbars
autoscrollbar(superscroll_obj, p);

end
