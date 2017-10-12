function h = plot_identity_line( ax )
%PLOT_IDENTITY_LINE 
% Plot identity line: y = x

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017
if nargin == 0
    ax = gca;
end
hold(ax,'on');

m = max(ax.XLim(1),ax.YLim(1));
M = min(ax.XLim(2),ax.YLim(2));
h = plot(ax, [m M], [m M], 'k--');
end

