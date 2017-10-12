function resize_max_square( ax, min_v )
%RESIZE_MAX_SQARE
% resize xlim and ylim to the largest containing interval

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

if nargin == 0
    ax = gca;
end

x = ax.XLim;
y = ax.YLim;

l = [min(x(1),y(1)), max(x(2),y(2))];

if nargin>0
    l(1) = min_v;
end

ax.XLim = l;
ax.YLim = l;
end

