function [labels, theta, phi] = read_chanlocs(chanlocfile, locExt)
% Robust parsing for common channel location file types.
labels = {};
theta = [];
phi = [];

fid = fopen(chanlocfile,'r');
if fid == -1
    error('Unable to open chanlocfile: %s', chanlocfile);
end
cleanup = onCleanup(@() fclose(fid));

switch locExt
    case 'xyz'
        % Expect: label X Y Z (cartesian)
        C = textscan(fid, '%s%f%f%f%*[^\n]', 'CommentStyle', 'c++', 'Delimiter', ' \t');
        if isempty(C{1}), error('No data read from .xyz file'); end
        labels = C{1};
        X = C{2}; Y = C{3}; Z = C{4};
        % cart2sph returns [azimuth elevation r] with azimuth in radians
        [az, el, ~] = cart2sph(X, Y, Z);
        theta = az * 180 / pi;   % azimuth -> theta (degrees)
        phi   = el * 180 / pi;   % elevation -> phi (degrees)
    case 'ced'
        % Common .ced format: index label theta radius X Y Z sph_theta sph_phi sph_radius
        % Skip header line if present
        frewind(fid);
        header = fgetl(fid);
        % Try to detect header by presence of non-numeric tokens
        if ischar(header) && ~isempty(regexp(header,'[A-Za-z]','once'))
            % header present, use textscan with headerlines=1
            frewind(fid);
            C = textscan(fid, '%n%s%f%f%f%f%f%f%f%f%*[^\n]', 'HeaderLines',1);
        else
            frewind(fid);
            C = textscan(fid, '%n%s%f%f%f%f%f%f%f%f%*[^\n]');
        end
        if isempty(C{2}), error('No data read from .ced file'); end
        labels = C{2};
        sph_theta = C{8};
        sph_phi   = C{9};
        % Adjust orientation consistent with original code
        theta = sph_theta + 90;
        theta(theta > 180) = theta(theta > 180) - 360;
        phi = sph_phi;
    case 'locs'
        % Expect: index theta radius label
        frewind(fid);
        C = textscan(fid, '%n%f%f%s%*[^\n]');
        if isempty(C{4}), error('No data read from .locs file'); end
        labels = C{4};
        th = C{2};
        radius = C{3};
        theta = -th + 90;
        theta(theta > 180) = theta(theta > 180) - 360;
        phi = 90 - (radius * 180);
    case 'csd'
        % Expect: label theta phi
        frewind(fid);
        C = textscan(fid, '%s%f%f%*[^\n]', 'CommentStyle', 'c++');
        if isempty(C{1}), error('No data read from .csd file'); end
        labels = C{1};
        theta = C{2};
        phi = C{3};
    otherwise
        error('Your channel location file extension ''.%s'' is not supported.', locExt);
end
end

