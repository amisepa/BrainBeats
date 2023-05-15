% Get channel neighors matrix from EEG channel locations
% using 2D polar (angular) projection and then 2D Delaunay triangulation.
% 
% inputs:
%       chanlocs - structure with EEG channel XYZ coordinates
%       vis - plot visualization (1) or not (0)
% 
% Cedric Cannard, Sep 2022

function [neighbors, channeighbstructmat] = get_channelneighbors(chanlocs,vis)

tmppath = fileparts(which('compute_mcc.m'));
addpath(fullfile(tmppath, 'subfunctions'))

cfg.elec.elecpos(:,1) = [ chanlocs.X ];
cfg.elec.elecpos(:,2) = [ chanlocs.Y ];
cfg.elec.elecpos(:,3) = [ chanlocs.Z ];
cfg.elec.label = { chanlocs.labels };
cfg.label = { chanlocs.labels };
cfg.method = 'triangulation';

% find channel neighbors
data = cfg;
cfg  = rmfield(cfg, 'label');     % first input must not be data
data = rmfield(data, 'method');   % second input must not be method

% neighbors = ft_prepare_neighbours(cfg, data); %%% TO EXTRACT
% cfg.feedback = 'no';
cfg.channel  = 'all';
cfg.compress = 'yes';
cfg.parcellation = 'parcellation';

% check if the input data is valid for this function
data = ft_checkdata(data);

% set the default for senstype depending on the data
cfg.senstype = ft_getopt(cfg, 'senstype', 'eeg');

% get 3D positions from the sensor description
sens = ft_fetch_sens(cfg, data);
chanpos = sens.chanpos;
label   = sens.label;

% remove channels that are not in data
[dum, sensidx] = match_str(data.label, label);
chanpos = chanpos(sensidx, :);
label   = label(sensidx);

% select the desired channels
desired = ft_channelselection(cfg.channel, label);
[sensidx] = match_str(label, desired);
chanpos = chanpos(sensidx, :);
label   = label(sensidx);

% Project sensor positions on 2D plane if not already the case
if size(chanpos, 2) == 2 || all( chanpos(:,3) == 0 )
    proj = chanpos(:,1:2); % already on a 2D plane
else
    % project sensor on a 2D plane (from function elproj)
    x = chanpos(:,1);
    y = chanpos(:,2);
    if size(chanpos, 2) == 3
        z = chanpos(:,3);
    end

    % use default polar (angular) projection
    [az, el, r] = cart2sph(x, y, z);
    [x, y] = pol2cart(az, pi/2 - el);
    proj = [x, y];
end

% 2D Delaunay triangulation of the projected points
tri = delaunay(proj(:,1), proj(:,2));
if strcmp(cfg.compress, 'yes')
    tri_x = delaunay(proj(:,1)./2, proj(:,2)); % compress in the x-direction
    tri_y = delaunay(proj(:,1), proj(:,2)./2); % compress in the y-direction
    tri = [tri; tri_x; tri_y];
end

% Compute the neighbourhood geometry from the gradiometer/electrode positions
% neighbors = compneighbstructfromtri(chanpos, label, tri);

% mark neighbors according to triangulation
nchan = length(label);
channeighbstructmat = zeros(nchan,nchan);
for i = 1:size(tri, 1)
    channeighbstructmat(tri(i, 1), tri(i, 2)) = 1;
    channeighbstructmat(tri(i, 1), tri(i, 3)) = 1;
    channeighbstructmat(tri(i, 2), tri(i, 1)) = 1;
    channeighbstructmat(tri(i, 3), tri(i, 1)) = 1;
    channeighbstructmat(tri(i, 2), tri(i, 3)) = 1;
    channeighbstructmat(tri(i, 3), tri(i, 2)) = 1;
end

% construct a structured cell-array with all neighbors
neighbors = struct;
alldist = [];
for i = 1:nchan
    neighbors(i).label          = label{i};
    neighbidx                   = find(channeighbstructmat(i,:));
    neighbors(i).dist           = sqrt(sum((repmat(chanpos(i, :), numel(neighbidx), 1) - chanpos(neighbidx, :)).^2, 2));
    alldist                     = [alldist; neighbors(i).dist];
    neighbors(i).neighblabel    = label(neighbidx);
end

% Remove neighbouring channels that are too far away (IMPORTANT if missing sensors)
maxdist = mean(alldist)+3*std(alldist);
for i = 1:nchan
    idx = neighbors(i).dist > maxdist;
    neighbors(i).dist(idx)         = [];
    neighbors(i).neighblabel(idx)  = [];
end
% neighbors = rmfield(neighbors, 'dist');

% Only select channels that are in the data
% if isfield(cfg, 'channel') && ~isempty(cfg.channel)
% %     desired = ft_channelselection(cfg.channel, data.label);
%     desired = data.label;
% end
% complete = struct;
% for i = 1:numel(desired)
% complete(i).label = desired{i};
% sel = find(strcmp({neighbors(:).label}, desired{i}));
% if numel(sel)==1
%   % take the set of neighbors from the definition
%   complete(i).neighblabel = neighbors(sel).neighblabel;
% else
%   % there are no neighbors defined for this channel
%   complete(i).neighblabel = {};
% end
% end
% neighbors = complete;

% Convert neighbors into row-arrays for a nicer code representation
% for i = 1:length(neighbors)
%   neighbors(i).neighblabel = neighbors(i).neighblabel(:)';
% end

% check that all chans have neighbors
k = 0;
for i = 1:length(neighbors)
    if isempty(neighbors(i).neighblabel)
        warning('no neighbors found for %s', neighbors(i).label);
    end
    k = k + length(neighbors(i).neighblabel);
end
if k==0
    warning('No neighbouring channels were specified or found');
else
    fprintf('there are on average %.1f neighbors per channel\n', k/length(neighbors));
end

% Visual feedback (or try convert_3dto2d for eeglab topo with nose and labels)
if vis
    tmpcfg = keepfields(cfg, {'layout', 'rows', 'columns', 'commentpos', 'skipcomnt', ...
        'scalepos', 'skipscale', 'projection', 'viewpoint', 'rotate', 'width', 'height', ...
        'elec', 'grad', 'opto', 'showcallinfo', 'trackcallinfo', 'trackconfig', ...
        'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
    tmpcfg.neighbours = neighbors;
    tmpcfg.senstype = cfg.senstype;
    ft_neighbourplot(tmpcfg, data);
end
