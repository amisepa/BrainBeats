function plot_topo(data,chanlocs,mode,dataType)
% PLOT_TOPO - Scalp topography of one value (or one curve) per EEG channel.
%
% Usage:
%   plot_topo(data, chanlocs, mode, dataType)
%
% Inputs:
%   data     - nChan x 1 values (mode 1), or nChan x nValues (mode 2, e.g.
%              multiscale entropy)
%   chanlocs - EEGLAB channel locations of the same nChan channels
%   mode     - 1 = 2D topoplot in the current axes (parula colormap, colorbar)
%              2 = 3D electrode plot in a new figure; with several values per
%                  channel, clicking an electrode plots its curve
%   dataType - 'psd' or 'entropy'. With 'entropy', values < .01 are set to
%              1e-4 and flagged as probable bad channels (mode 1).
%
% Copyright (C) - Cedric Cannard, 2023

chanlabels = {chanlocs.labels};
x = [ chanlocs.X ]';
y = [ chanlocs.Y ]';
z = [ chanlocs.Z ]';

% Rotate X Y Z coordinates (for the 3D plot; topoplot uses chanlocs directly)
% rotate = 0;       %nosedir = +x
rotate = 3*pi/2;    %nosedir = +y
% rotate = pi;      %nosedir = -x
% rotate = pi/2;
allcoords = (y + x.*sqrt(-1)).*exp(sqrt(-1).*rotate);
x = imag(allcoords);
y = real(allcoords);

% Project 3D positions on 2D plane if not already done
chanpos(:,1) = x;
chanpos(:,2) = y;
chanpos(:,3) = z;

if all(chanpos(:,3)==0)
    coord = chanpos(:,1:2); % positions already projected on a 2D plane
else
    coord = chanpos; % use 3-D data for plotting
end

%% 2D scalp topography with colors (band-power, IAF, entropy)
if mode == 1

    % Deal with near-0 values (for entropy data)
    if strcmp(dataType,'entropy')
        idx = data < 0.01;
        if sum(idx) > 0
            data(idx) = 0.0001;
            warning(['Channel ' chanlocs(idx).labels ' is probably a bad channel.'])
        end
    end

    % Scalp topo (color limits = data range; fails silently if all values are equal)
    topoplot(data, chanlocs,'emarker',{'.','k',7,1},'electrodes','on');
    try
        set(gca,'CLim',[min(data) max(data)]);   % clim needs R2022a
    catch
    end
    colormap('parula');
    c = colorbar;
    c.Label.FontSize = 12;
    c.Label.FontWeight = 'bold';
end

%% 3D electrode plot; click an electrode to plot its values (nChan x nValues data)
if mode == 2

    p = figure('color','w');
    axis equal
    axis vis3d
    axis off
    hold on

    for iChan = 1:size(data,1)

        if length(data(iChan,:)) == 1 % measures with one value per channel
            % electrode marker
            p(iChan) = plot3(coord(iChan,1),coord(iChan,2),coord(iChan,3), ...
                'MarkerEdgeColor','k','MarkerFaceColor', 'k', ...
                'Marker','o','MarkerSize',5);

            % Display channel label + value for each channel
            text(coord(iChan,1)-15,coord(iChan,2)+10,coord(iChan,3), ...
                sprintf('%s: %6.1f',chanlabels{iChan}, ...
                round(data(iChan,:),2)),'FontSize',10,'fontweight','bold');

        else % several values per channel (e.g., multiscale entropy): clickable marker
            p(iChan) = plot3(coord(iChan,1),coord(iChan,2),coord(iChan,3), ...
                'MarkerEdgeColor','k','MarkerFaceColor', 'k', ...
                'Marker','o','MarkerSize', 5, 'UserData',iChan, ...
                'ButtonDownFcn', @(~,~,~) buttonCallback(data(iChan,:), coord(iChan,:), chanlabels{iChan}));

            % Display channel label above each electrode
            text(coord(iChan,1)-7,coord(iChan,2)+10,coord(iChan,3), ...
                sprintf('%s %6.3f',chanlabels{iChan}), ...
                'FontSize',10,'fontweight','bold');
            title('[Click on sensors to display entropy values]', ...
                'Position', [1 120 1], 'fontweight', 'bold')

        end
    end
end



%% subfunction to display data on click
function buttonCallback(tmpdata, coor, label)

% Plot the clicked channel's values across time scales
figure('color','w','Position', [500 500 280 210]);
plot(tmpdata,'linewidth',2,'color','black');
title(label,'FontSize',14)
xlabel('Time scale','FontSize',12,'fontweight','bold');
ylabel('Entropy','FontSize',12,'fontweight','bold')
