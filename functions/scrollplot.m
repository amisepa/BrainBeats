function scrollplot(varargin)
    % Extract inputs
    toPlot1 = varargin{1};      % Data to plot
    scrollOrient = varargin{2}; % 'X' for horizontal and 'Y' for vertical
    toDisplay = varargin{3};    % Window size (in seconds)

    toPlot2 = [];
    toPlot3 = [];
    if nargin >= 4, toPlot2 = varargin{4}; end
    if nargin == 5, toPlot3 = varargin{5}; end

    % Convert window size to percentage
    toDisplay = toDisplay / max(toPlot1{1});

    % Plot main data
    p = plot(toPlot1{:}, 'linewidth', 1);
    ax = p.Parent;
    XX = p.XData;
    YY = p.YData;

    % Compute x-axis limits
    totxlim = [min(XX, [], 'omitnan'), max(XX, [], 'omitnan')];
    xlim_init = [totxlim(1), totxlim(1) + toDisplay * diff(totxlim)];
    set(ax, 'XLim', xlim_init);
    updateYLim(ax, XX, YY);  % Set initial Y-limits

    hold on
    if ~isempty(toPlot2), plot(toPlot2{:}); end
    if ~isempty(toPlot3), plot(toPlot3{:}); end

    % Scroll initialization
    superscroll_obj = superscroll(ax, scrollOrient);
    autoscrollbar(superscroll_obj, p);

    % Add listener to clamp XLim and update YLim
    addlistener(ax, 'XLim', 'PostSet', @(src, evt) handleScroll(ax, XX, YY, totxlim, toDisplay));

    % === Nested Functions ===
    function updateYLim(ax, XData, YData)
        xlims = get(ax, 'XLim');
        visible = XData >= xlims(1) & XData <= xlims(2);
        Yvis = YData(visible);
        if ~isempty(Yvis) && any(~isnan(Yvis))
            margin = 0.1 * range(Yvis);
            if margin == 0, margin = 1; end
            set(ax, 'YLim', [min(Yvis)-margin, max(Yvis)+margin]);
        end
    end

    function handleScroll(ax, XData, YData, totxlim, toDisplayFrac)
        winWidth = toDisplayFrac * diff(totxlim);
        cur_xlim = get(ax, 'XLim');
        % Clamp
        if cur_xlim(1) < totxlim(1)
            cur_xlim = [totxlim(1), totxlim(1) + winWidth];
        elseif cur_xlim(2) > totxlim(2)
            cur_xlim = [totxlim(2) - winWidth, totxlim(2)];
        end
        set(ax, 'XLim', cur_xlim);
        updateYLim(ax, XData, YData);
    end
end
