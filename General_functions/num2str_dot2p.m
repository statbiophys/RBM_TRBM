function [ s ] = num2str_dot2p( x )
%NUM2STR_DOT2P as num2str, but replaces the dot with a p

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

s = num2str(x);
s(s=='.') = 'p';
s(s=='-') = 'm';
end

