function [ Pv_l ] = RBM_h2Pv( M, h )
%RBM_V2PH
% M : model
% h : vector column of hidden states
% v : vector column of mean visible units

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

Pv_l = sigmoid_approx( bsxfun(@plus, M.w'*h, M.a') );

end

