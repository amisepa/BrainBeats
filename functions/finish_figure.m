% FINISH_FIGURE - Final touches on BrainBeats figures, and draw them now.
%
% Hides the figure toolbar and the axes toolbars (shown on every axes in
% recent MATLAB versions), undocks a figure that opened docked when the
% default window style is not docked, and draws the figure immediately:
% since R2025a figures are only rendered when MATLAB flushes the graphics
% queue, so without this they stay blank during the computations that follow
% (e.g. ICA).
%
% Usage:
%   finish_figure            % current figure
%   finish_figure(figs)      % figure handle(s)
%
% Copyright (C) - Cedric Cannard, 2026

function finish_figure(figs)

if nargin < 1, figs = get(groot,'CurrentFigure'); end
figs = figs(isgraphics(figs,'figure'));
for fig = figs(:)'
    try
        set(fig, 'ToolBar', 'none');
        for ax = findall(fig, 'Type', 'axes')'
            if isprop(ax,'Toolbar') && ~isempty(ax.Toolbar)
                ax.Toolbar.Visible = 'off';
            end
        end
        if strcmp(get(fig,'WindowStyle'),'docked') && ~strcmp(get(groot,'DefaultFigureWindowStyle'),'docked')
            set(fig, 'WindowStyle', 'normal');
        end
    catch
        % older MATLAB versions (no axes toolbar) or figures that refuse the change
    end
end
drawnow
