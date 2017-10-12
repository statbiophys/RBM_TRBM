function [ Y ] = circshift_asOLD( X, K, D )
%CIRCSHIFT_ASOLD
% for circshift as in old versions of Matlab

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

if D == 1
    Y = circshift(X,K);
else
    o_l = 1:length(size(X));
    o_l([1 D]) = [D 1];
    Y = permute(circshift(permute(X,o_l), K), o_l);
end

% if 0
%     %% test
%     X = rand(5,8,6);
%     K = 3;
%     D = 2;
%     Y1 = circshift(X, K, D);
%     Y2 = circshift_asOLD(X, K, D);
%     isequal(Y1, Y2)
% end
end

