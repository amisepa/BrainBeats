function scrollplot(varargin)
    % scrollplot - Scrollable time series plot with optional overlays.
    %
    % Scrolling via:
    %   - Slider bar at the bottom of the axes
    %   - Left/Right arrow keys (pan by 50% of window)
    %   - Shift + arrow keys (pan by 10% of window, fine control)
    %
    % Usage:
    %   scrollplot(toPlot1, scrollOrient, windowSize)
    %   scrollplot(toPlot1, scrollOrient, windowSize, toPlot2, toPlot3, ...)
    %
    % Inputs:
    %   toPlot1      - {time, signal, 'prop', val, ...} for main signal
    %   scrollOrient - {'X'} (only horizontal supported)
    %   windowSize   - visible window width in x-axis units (e.g. seconds)
    %   toPlot2, ... - (optional) additional overlay series as cell arrays

    toPlot1      = varargin{1};
    % scrollOrient = varargin{2};  % always X, kept for API compatibility
    win_size     = varargin{3};
    extraPlots   = {};
    if nargin > 3
        extraPlots = varargin(4:end);
    end

    % Enforce column vectors
    x = toPlot1{1}(:);
    y = toPlot1{2}(:);
    if length(x) ~= length(y)
        error('scrollplot: X and Y must have the same number of elements.');
    end
    toPlot1{1} = x;
    toPlot1{2} = y;

    % Axis limits
    x_min    = min(x, [], 'omitnan');
    x_max    = max(x, [], 'omitnan');
    win_size = min(win_size, x_max - x_min);  % clamp to data range

    % Plot main signal
    p  = plot(toPlot1{:}, 'LineWidth', 1);
    ax = p.Parent;
    fig = ax.Parent;

    set(ax, 'XLim', [x_min, x_min + win_size]);
    updateYLim();
    hold on

    % Plot overlays
    for k = 1:length(extraPlots)
        if ~isempty(extraPlots{k})
            plot(extraPlots{k}{:});
        end
    end

    % --- Slider ---
    % Placed in normalized figure coords below the current axes position
    ax_pos = get(ax, 'Position');  % [left bottom width height] normalized
    slider_h = 0.03;
    slider = uicontrol(fig, 'Style', 'slider', ...
        'Units',    'normalized', ...
        'Position', [ax_pos(1), ax_pos(2) - slider_h - 0.005, ax_pos(3), slider_h], ...
        'Min',      x_min, ...
        'Max',      max(x_min + 1e-6, x_max - win_size), ...
        'Value',    x_min, ...
        'Callback', @sliderCallback);

    % Shrink axes slightly to make room for slider
    set(ax, 'Position', [ax_pos(1), ax_pos(2) + slider_h, ax_pos(3), ax_pos(4) - slider_h]);

    % --- Keyboard scrolling ---
    set(fig, 'KeyPressFcn', @keyCallback);

    % === Callbacks ===
    function sliderCallback(src, ~)
        t0 = src.Value;
        set(ax, 'XLim', [t0, t0 + win_size]);
        updateYLim();
    end

    function keyCallback(~, evt)
        step_large = win_size * 0.9;   % arrow key: 90% step (10% overlap), e.g. 1-10s -> 9-18s
        step_small = win_size * 0.1;   % shift+arrow: 10% step (fine control)
        cur = get(ax, 'XLim');
        t0  = cur(1);

        if strcmp(evt.Modifier, 'shift')
            step = step_small;
        else
            step = step_large;
        end

        switch evt.Key
            case 'rightarrow'
                t0 = min(t0 + step, x_max - win_size);
            case 'leftarrow'
                t0 = max(t0 - step, x_min);
            otherwise
                return
        end

        set(ax, 'XLim', [t0, t0 + win_size]);
        set(slider, 'Value', max(slider.Min, min(slider.Max, t0)));
        updateYLim();
    end

    function updateYLim()
        xlims   = get(ax, 'XLim');
        visible = x >= xlims(1) & x <= xlims(2);
        Yvis    = y(visible);
        if ~isempty(Yvis) && any(~isnan(Yvis))
            ymin = min(Yvis);
            ymax = max(Yvis);
            pp   = ymax - ymin;
            % if pp > 200
                set(ax, 'YLim', [ymin - 0.2*pp, ymax + 0.2*pp]);
            % else
            %     set(ax, 'YLim', [-150, 100]);
            % end
        end
    end
end