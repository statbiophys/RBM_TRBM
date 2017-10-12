function [ l ] = plot_hline( y, col, ax, line_style )
%PLOT_VLINE
% Plot vertical line

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017
if nargin<2
    col = 'k';
end
if nargin<3
    ax = gca;
end
if nargin<4
    line_style = '-';
end
hold(ax,'on');
l = plot(ax.XLim, [y y],'Color',col,'LineStyle', line_style);
end

