function [ Pv_l, mE_l] = TRBM_cyclic_h2Pv( M, h )
%RBM_V2PH
% M : model
% h : vector column of hidden states
% Pv_l : vector column of mean visible units

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

Pv_l = TRBM_cyclic_h2mEv( M, h );
if nargout>1
    mE_l = Pv_l; % mE of the +1 state
end
Pv_l = sigmoid_approx( Pv_l );

end

