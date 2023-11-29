% Get channel neighors matrix from EEG channel locations
% using 2D polar (angular) projection and then 2D Delaunay triangulation.
% 
% inputs:
%       chanlocs - structure with EEG channel XYZ coordinates
%       compress - compresssion (1, more neighbors, default) or not (0, less neighbors)
% 
% Example:
%       [neighbors, channeighbstructmat] = get_channelneighbors(chanlocs,compress)
% 
% Cedric Cannard, Sep 2022

function [neighbors, channeighbstructmat] = get_channelneighbors(chanlocs)

compress = true; 

% Project sensor positions on 2D plane
x = [chanlocs.X]';
y = [chanlocs.Y]';
z = [chanlocs.Z]';

% use default polar (angular) projection
[az, el, r] = cart2sph(x, y, z);
[x, y] = pol2cart(az, pi/2 - el);
proj = [x, y];

% 2D Delaunay triangulation of the projected points
tri = delaunay(proj(:,1), proj(:,2));
if compress
    tri_x = delaunay(proj(:,1)./2, proj(:,2)); % compress in the x-direction
    tri_y = delaunay(proj(:,1), proj(:,2)./2); % compress in the y-direction
    tri = [tri; tri_x; tri_y];
end

% mark neighbors according to triangulation
nchan = length(x);
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
chanpos = [x(:) y(:)];
neighbors = struct;
alldist = [];
for i = 1:nchan
    neighbors(i).label          = chanlocs(i).labels;
    neighbidx                   = find(channeighbstructmat(i,:));
    neighbors(i).dist           = sqrt(sum((repmat(chanpos(i, :), numel(neighbidx), 1) - chanpos(neighbidx, :)).^2, 2));
    alldist                     = [alldist; neighbors(i).dist];
    neighbors(i).neighblabel    = {chanlocs(neighbidx).labels};
end

% Remove neighbouring channels that are too far away (IMPORTANT if missing sensors)
maxdist = mean(alldist)+3*std(alldist);
for i = 1:nchan
    idx = neighbors(i).dist > maxdist;
    neighbors(i).dist(idx)         = [];
    neighbors(i).neighblabel(idx)  = [];
end
% neighbors = rmfield(neighbors, 'dist');

% Convert neighbors into row-arrays for a nicer code representation
for i = 1:length(neighbors)
  neighbors(i).neighblabel = neighbors(i).neighblabel(:)';
end

% check that all chans have neighbors and report average
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
    fprintf('There are on average %.1f neighbors per channel\n', k/length(neighbors));
end

% Visual feedback (or try convert_3dto2d for eeglab topo with nose and labels)
% if vis
%     tmpcfg = keepfields(cfg, {'layout', 'rows', 'columns', 'commentpos', 'skipcomnt', ...
%         'scalepos', 'skipscale', 'projection', 'viewpoint', 'rotate', 'width', 'height', ...
%         'elec', 'grad', 'opto', 'showcallinfo', 'trackcallinfo', 'trackconfig', ...
%         'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
%     tmpcfg.neighbours = neighbors;
%     tmpcfg.senstype = cfg.senstype;
%     ft_neighbourplot(tmpcfg, data);
% end
