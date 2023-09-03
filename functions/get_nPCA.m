%% Get the number of principal components that should be retained to
% explain X% of variance in the data. This is useful for
% diemnsion-reduction to use the minimum number of components to retain. 
% 
% Inputs:
%   data
%   varThresh   - variance threshold to retain in % (default = 99)
%   vis         - visualize (true) or not (false)
% 
% Usage: 
%   nComps = get_nPCA(data, varThresh, vis)
%   nComps = get_nPCA(data, 99, true)
% 
% Copyright (C) - Cedric Cannard, 2023

function nComps = get_nPCA(data, varThresh, vis)

if ~exist('varThresh','var') || isempty(varThresh)
    varThresh = 99;
end
if ~exist('varThresh','var') || isempty(vis)
    vis = 1;
end

% Calculate cumulative variance explained
[coeff, score, ~, ~, explained] = pca(double(data)); 
cumulativeExplained = cumsum(explained);

% number of components explaining 99% of variance
nComps = find(cumulativeExplained >= varThresh, 1); 

% Plot cumulative variance explained
if vis
    figure('color','w');
    plot(cumulativeExplained, 'o-');
    xlabel('Number of Components');
    ylabel('Cumulative Variance Explained (%)');
    title(sprintf('# of components to retain to preserve 99%% of variance = %g', nComps));
    grid on; ylim([0 100])
end


%% some other PCA stuff

% X = table2array(HRV);
% % X = zscore(X); % normalize
% 
% [coeff,score,latent,~,explained] = pca(X); % uses Singular Value Decomposition (default)
% 
% % % Calculate eigenvalues and eigenvectors of the covariance matrix
% % covarianceMatrix = cov(X);
% % [V,D] = eig(covarianceMatrix);
% % coeff
% % V
% 
% % score shows projections of the original data on the principal component vector space.
% %  (same as doing X*coeff)
% 
% % The columns of score are orthogonal to each other.
% % corrcoef(score)
% 
% % The variances of these vectors are the eigenvalues of the covariance matrix,
% % and are also the output "latent".
% % var(score)'
% % latent
% 
% % Now you can think about dimension reduction. Take a look at the variable 
% % 'explained'. It tells you how much of the variation is captured by each column 
% % of 'score'. Here is where you have to make a judgement call. How much 
% % of the total variation are you willing to ignore? 
% % One guideline is that if you plot explained, there will often be an "elbow" 
% % in the plot, where each additional variable explains very little additional 
% % variation. Keep only the components that add a lot more explanatory power, 
% % and ignore the rest.
% % If first 3 components together explain 87% of the variation; suppose 
% % that's good enough. Then, you would only keep those 3 dimensions (the first 3 
% % columns of score). You will have 7 observations in 3 dimensions (variables) instead of 5.
% 
% % explained' shows the variation catpured by each column in 'score'
% 
% % make a biplot (PC1 on x-axis and PC2 on y-axis)
% xPC = 1;
% yPC = 2;
% figure; hold on
% for nc = 1:N
%     plot([0 coeff(nc,xPC)],[0 coeff(nc,yPC)],'b'); 
% end
% 
% % Scatter plot of the first two components (normalize variables original data by dividng by SD)
% % figure; hold on
% % h = plot(score(:,1),score(:,2),'ro');
% % set(h,'MarkerSize',16)
% % set(gca,'XLim',[-1 1],'YLim',[-1 1],'Box','on')
% % axis square
% % xlabel('Component 1')
% % ylabel('Component 2')
% % % Add a circle
% % p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 0.8);
% % plot(p, 'FaceColor', 'w', 'EdgeColor', 'r')
% 
% % Recover the original data from the PC
% % score * coeff'
% 
% % Default (normalized varimax) rotation: first 3 principal components.
% % [1] Harman, H. H. Modern Factor Analysis. 3rd ed. Chicago: University of Chicago Press, 1976.
% % [2] Lawley, D. N., and A. E. Maxwell. Factor Analysis as a Statistical Method. 2nd ed. New York: American Elsevier Publishing, 1971. 
% % LPC = pca(X);
% [L1,T] = rotatefactors(coeff(:,1:2));
% inv(T'*T) % Correlation matrix of the rotated factors
