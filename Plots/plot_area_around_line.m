function h = plot_area_around_line( x_l, y_l, e_l , col, alpha, Ymin)
%PLOT_AREA_AROUND_CURVE

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

lower_l = y_l(:)-e_l(:);
upper_l = flipud(y_l(:)+e_l(:));
if nargin>5
    lower_l = max(lower_l, Ymin);
    upper_l = max(upper_l, Ymin);
end

h = patch([x_l(:); flipud(x_l(:))],[lower_l; upper_l],col);

if 1
    set(h,'edgecolor','none');
else
    set(h,'edgecolor',col);
end

if nargin>4 && alpha ~= 1
    set(h,'facealpha',alpha)
end
end

