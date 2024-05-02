function plot_topo(data,chanlocs,mode,dataType)

chanlabels = {chanlocs.labels};
x = [ chanlocs.X ]';
y = [ chanlocs.Y ]';
z = [ chanlocs.Z ]';

% Rotate X Y Z coordinates
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

    % Scalp topo
    % if strcmpi(dataType,'entropy') && var(data) < 0.1
    %     disp('Entropy data do not have enough variance across electrodes to plot an informative scalp topography.')
    % else
    % figure('color','w');
    topoplot(data, chanlocs,'emarker',{'.','k',7,1},'electrodes','on');
    try
        clim([min(data) max(data)]);
    catch
    end
    % if strcmpi(dataType,'psd')
    colormap('parula');
    % else
    %     colormap('hot');
    % end
    c = colorbar;
    % ylabel(c,'Coherence','FontSize',12)
    % title('Entropy','FontSize',10);
    % if strcmpi(dataType,'psd')
    % c.Label.String = 'Power';
    % elseif strcmpi(dataType,'iaf')
    %     c.Label.String = 'Alpha centre of gravity';
    % elseif strcmpi(dataType,'entropy')
    %     c.Label.String = 'Entropy';
    % end
    c.Label.FontSize = 12;
    c.Label.FontWeight = 'bold';
    % end
end

%% 3D scalp topo allowing to open data for each electrode (for 2D data)
if mode == 2

    p = figure('color','w');
    % p.Position = [100 100 540 400];
    axis equal
    axis vis3d
    axis off
    hold on

    % adj = mean(data(1,:))*5; % to scale marker size

    for iChan = 1:size(data,1)

        if length(data(iChan,:)) == 1 % measures with one value per channel
            % 3D plot of entropy values at electrode locations
            p(iChan) = plot3(coord(iChan,1),coord(iChan,2),coord(iChan,3), ...
                'MarkerEdgeColor','k','MarkerFaceColor', 'k', ...
                'Marker','o','MarkerSize',5);

            % Display channel label + entropy value for each channel
            text(coord(iChan,1)-15,coord(iChan,2)+10,coord(iChan,3), ...
                sprintf('%s: %6.1f',chanlabels{iChan}, ...
                round(data(iChan,:),2)),'FontSize',10,'fontweight','bold');

        else % for multiscales, take area under the curve as sensor size
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

% Entropy measures with only one value per channel
figure('color','w','Position', [500 500 280 210]);
plot(tmpdata,'linewidth',2,'color','black'); % blue: [0, 0.4470, 0.7410]
% area(tmpdata,'linewidth',2);
title(label,'FontSize',14)
% xticks(2:nScales); xticklabels(join(string(scales(:,2:end)),1)); xtickangle(45)
% xlim([2 nScales]);
xlabel('Time scale','FontSize',12,'fontweight','bold');
ylabel('Entropy','FontSize',12,'fontweight','bold')
