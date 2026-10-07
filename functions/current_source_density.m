function [X, Y] = current_source_density(Data, G, H, lambda, head)
% Vectorized implementation of the CSD algorithm (Perrin / Kayser).
% Data: nElec x nSamples
% G, H: nElec x nElec
% lambda: smoothing constant
% head: head radius in cm (if provided), scaling uses head^2 as in original code
% 
% CSD - Current Source Density (CSD) transformation based on spherical spline
%       surface Laplacian as suggested by Perrin et al. (1989, 1990)
%
% (published in appendix of Kayser J, Tenke CE, Clin Neurophysiol 2006;117(2):348-368)
%
% Usage: [X, Y] = CSD(Data, G, H, lambda, head);
%
% Implementation of algorithms described by Perrin, Pernier, Bertrand, and
% Echallier in Electroenceph Clin Neurophysiol 1989;72(2):184-187, and 
% Corrigenda EEG 02274 in Electroenceph Clin Neurophysiol 1990;76:565.
% 
% Input parameters:
%   Data = surface potential electrodes-by-samples matrix
%      G = g-function electrodes-by-electrodes matrix
%      H = h-function electrodes-by-electrodes matrix
% lambda = smoothing constant lambda (default = 1.0e-5)
%   head = head radius (default = no value for unit sphere [µV/m²])
%          specify a value [cm] to rescale CSD data to smaller units [µV/cm²]
%          (e.g., use 10.0 to scale to more realistic head size)
%
% Output parameter:
%      X = current source density (CSD) transform electrodes-by-samples matrix
%      Y = spherical spline surface potential (SP) interpolation electrodes-
%          by-samples matrix (only if requested)
%
% Copyright (C) 2003 by Jürgen Kayser (Email: kayserj@pi.cpmc.columbia.edu)
% GNU General Public License (http://www.gnu.org/licenses/gpl.txt)
% Updated: $Date: 2005/02/11 14:00:00 $ $Author: jk $
%        - code compression and comments 
% Updated: $Date: 2007/02/07 11:30:00 $ $Author: jk $
%        - recommented rescaling (unit sphere [µV/m²] to realistic head size [µV/cm²])
%   Fixed: $Date: 2009/05/16 11:55:00 $ $Author: jk $
%        - memory claim for output matrices used inappropriate G and H dimensions
%   Added: $Date: 2009/05/21 10:52:00 $ $Author: jk $
%        - error checking of input matrix dimensions
%

[nElec, nPnts] = size(Data);

% Basic dimension checks
if ~(size(G,1) == size(G,2)) || ~(size(H,1) == size(H,2)) || ...
   ~(size(G,1) == nElec) || ~(size(H,1) == nElec)
    X = NaN; Y = NaN;
    error('G and H must be nElec-by-nElec and match rows of Data.');
end

% Center data (remove grand mean across electrodes for each sample)
mu = mean(Data,1);            % 1 x nPnts
Z = Data - mu;                % nElec x nPnts

% Defaults
if nargin < 5 || isempty(head), head = 1.0; end
if nargin < 4 || isempty(lambda), lambda = 1.0e-5; end

% Add smoothing constant to diagonal of G
G = G + lambda * eye(nElec);

% Solve for Gi efficiently (we need Gi * Z). Use factorization/backslash.
% Compute Gi = inv(G) via backslash on identity (numerically stable)
Gi = G \ eye(nElec);         % nElec x nElec

% Precompute TC (row sums of Gi) and sgi
TC = sum(Gi,2);              % nElec x 1
sgi = sum(TC);               % scalar

% Compute Cp for all samples at once
Cp = Gi * Z;                 % nElec x nPnts

% Compute c0 for each sample (1 x nPnts)
c0 = sum(Cp,1) ./ sgi;       % 1 x nPnts

% Compute C matrix: C = Cp - TC * c0
C = Cp - TC * c0;            % nElec x nPnts

% Scale factor for head radius: original code squares head
headScale = head * head;
if headScale == 0, headScale = 1; end

% Compute X (CSD) and Y (spherical potential interpolation) vectorized
X = (H * C) ./ headScale;    % nElec x nPnts
Y = G * C + repmat(c0, nElec, 1); % nElec x nPnts

end
