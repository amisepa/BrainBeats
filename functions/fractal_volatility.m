function [dimension, standard_dev] = fractal_volatility(data)
% FRACTAL_VOLATILITY - Fractal dimension of a 1-D signal by box counting.
%
% The signal is treated as a function of its indices and rescaled into the
% unit square (x = sample index, y = amplitude). Boxes of width 2^-j are
% counted down to the sampling resolution, and the dimension is the OLS slope
% of log(count) vs. log(1/width), after discarding scales whose local slope
% deviates from the median by more than IQR/2. Adapted from fractalvol.
%
% Usage:
%   [dimension, standard_dev] = fractal_volatility(data)
%
% Inputs:
%   data         - signal (vector), or N x 2 matrix [x y]
%
% Outputs:
%   dimension    - fractal dimension, rounded to 3 decimals (1 = smooth line,
%                  up to 2 = plane-filling)
%   standard_dev - uncertainty of the slope, sqrt(SSE * inv(X'X))(2,2)

if size(data,1) < size(data,2)
    data = data';
end

if size(data,2) == 1
    data = [(1:length(data))' data];
end

max1 = max(data(:,1));
min1 = min(data(:,1));
max2 = max(data(:,2));
min2 = min(data(:,2));

%normalize all to unit square
normdata = data;
normdata(:,1) = normdata(:,1)-min1;
normdata(:,1) = normdata(:,1)./(max1-min1);
normdata(:,2) = normdata(:,2)-min2;
normdata(:,2) = normdata(:,2)./(max2-min2);

% smallest box width stays above the sample spacing, so that every box
% column contains at least one sample
minwidth = min(diff(normdata(:,1)));
minwidth = log2(minwidth);
minwidth = abs(ceil(minwidth));
minwidth = minwidth-1;

n = zeros(minwidth,1);

% box count for each width 2^-j
    for j = 1:minwidth
        width = 2^-j;
        xaxis_pos = 0;
        boxcount = 0;
        while xaxis_pos<1
            indx = (xaxis_pos <= normdata(:,1) &...
                normdata(:,1)<xaxis_pos+width);
            if 1-xaxis_pos == width
                indx(end) = true;
            end
            
            vertical_column = normdata(indx,2);
            
            if length(vertical_column) == 1
                boxcount = boxcount + 1;
            else
                rawcount = (max(vertical_column)-min(vertical_column))/width;
                rawcount = rawcount + rem(min(vertical_column),width);
                count = ceil(rawcount);
                boxcount = boxcount + count;
            end
            xaxis_pos = xaxis_pos + width; %advance on x axis
        end
        n(j) = boxcount;
    end

r = 2.^-(1:minwidth);
r = r';

% local slopes; discard scales that deviate from the median slope by more than IQR/2
s=-gradient(log(n))./gradient(log(r));
IQR = iqr(s);

indx2 = abs(s-median(s)) > IQR/2;

x2 = log(r);
y2 = log(n);

s(indx2)= [];
x2(indx2) = [];
y2(indx2) = [];

% OLS fit of log(count) on log(width): dimension = -slope, with its uncertainty

X = [ones(size(x2)) x2];
beta = pinv(X)*y2;
C = pinv(X'*X);
e=y2-X*pinv(X)*y2;
sigma = e'*e*C;
sigma = sqrt(sigma);

dimension = round(-beta(2),3);
standard_dev = sigma(2,2);
