function scrollplot(varargin)
    % Extract inputs
    toPlot1 = varargin{1};      % Data to plot: {time, signal, 'color', ...}
    scrollOrient = varargin{2}; % 'X' for horizontal or 'Y' for vertical scroll
    toDisplay = varargin{3};    % Window size (in seconds or units of X-axis)

    toPlot2 = [];
    toPlot3 = [];
    if nargin >= 4, toPlot2 = varargin{4}; end
    if nargin == 5, toPlot3 = varargin{5}; end

    % Enforce shape compatibility for primary data
    x = toPlot1{1}(:);  % force column vector
    y = toPlot1{2}(:);  % force column vector
    if length(x) ~= length(y)
        error('X and Y inputs must have the same number of elements.');
    end
    toPlot1{1} = x;
    toPlot1{2} = y;

    % Convert window size to fraction
    toDisplay = toDisplay / max(x);

    % Plot main data
    p = plot(toPlot1{:}, 'LineWidth', 1);
    ax = p.Parent;
    XX = p.XData;
    YY = p.YData;

    % Compute initial X-limits
    totxlim = [min(XX, [], 'omitnan'), max(XX, [], 'omitnan')];
    xlim_init = [totxlim(1), totxlim(1) + toDisplay * diff(totxlim)];
    set(ax, 'XLim', xlim_init);
    updateYLim(ax, XX, YY);

    hold on
    if ~isempty(toPlot2), plot(toPlot2{:}); end
    if ~isempty(toPlot3), plot(toPlot3{:}); end

    % Scroll initialization
    superscroll_obj = superscroll(ax, scrollOrient);
    autoscrollbar(superscroll_obj, p);

    % Add listener for axis updates
    addlistener(ax, 'XLim', 'PostSet', @(src, evt) handleScroll(ax, XX, YY, totxlim, toDisplay));

    % === Nested Functions ===
    function updateYLim(ax, XData, YData)
        xlims = get(ax, 'XLim');
        visible = XData >= xlims(1) & XData <= xlims(2);
        Yvis = YData(visible);

        if ~isempty(Yvis) && any(~isnan(Yvis))
            ymin = min(Yvis);
            ymax = max(Yvis);
            peak_to_peak = ymax - ymin;

            if peak_to_peak > 200  % autoscale
                margin = 0.1 * peak_to_peak;
                set(ax, 'YLim', [ymin - margin, ymax + margin]);
            else  % fixed range
                set(ax, 'YLim', [-100, 100]);
            end
        end
    end

    % function updateYLim(ax, XData, YData)
    %     if isempty(YData) || all(isnan(YData))
    %         return;
    %     end
    %     y_min = min(YData, [], 'omitnan');
    %     y_max = max(YData, [], 'omitnan');
    %     margin = 0.05 * (y_max - y_min);
    %     set(ax, 'YLim', [y_min - margin, y_max + margin]);
    % end

    function handleScroll(ax, XData, YData, totxlim, toDisplayFrac)
        winWidth = toDisplayFrac * diff(totxlim);
        cur_xlim = get(ax, 'XLim');
        % Clamp X-limits
        if cur_xlim(1) < totxlim(1)
            cur_xlim = [totxlim(1), totxlim(1) + winWidth];
        elseif cur_xlim(2) > totxlim(2)
            cur_xlim = [totxlim(2) - winWidth, totxlim(2)];
        end
        set(ax, 'XLim', cur_xlim);
        updateYLim(ax, XData, YData);

    end
end
