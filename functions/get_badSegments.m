% Get latency bounds of artifactual segments
% 
% Inputs:
%   - mask: array of zeros (good sample) and ones (bad samples)
%   - minSize: minimum size (in samples) of artifactiual segment to be 
%           considered artifact (default = 0.05*srate, i.e. 50 ms)
%   - minGap: when gap size (in samples) between segments is smaller, they
%           are merged into a larger epoch to be removed to avoid many
%           discontinuities (fefault = 0.1*srate, i.e. 100 ms)
 
% Output:
%   - out: bounds of each bad segment (in samples)
% 
% Cedric Cannard, 2022

function out = get_badSegments(mask, minSize, minGap)

% Segment bounds
out = reshape(find(diff([false mask false])),2,[])';
out(:,2) = out(:,2)-1;

% Remove segments smaller than art_size
tmpSize = out(:,2) - out(:,1);
out(tmpSize < minSize,:) = [];

% Find segments with gaps between them smaller than minGap
if size(out,1) > 1
    for iSeg = 2:length(out)
        gapSize(iSeg,:) = out(iSeg,1) - out(iSeg-1,2);
    end
    gapSize(1,:) = [];
    idx = gapSize < minGap;
    
    % Merge them
    tmpout = out;
    for iSeg = 2:length(out)-1
        if idx(iSeg)
            out(iSeg,2) = out(iSeg+1,2);
            out(iSeg+1,:) = NaN;
        end
    end
    out(isnan(out(:,1)),:) = [];
end
