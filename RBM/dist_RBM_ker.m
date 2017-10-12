function [ d ] = dist_RBM_ker( x, y, ker, pow )
%DIST_RBM_KER

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

% d = x(:) - y(:);
d = x - y;
if size(d,1) ~= size(ker,2)
    d = d(:);
end
d = d.*(ker*d);
if nargin<4
    %     pow = 1;
    d = sum(d(:));
else
    d = sum(d(:).^pow);
end

d = sqrt(d);
end

