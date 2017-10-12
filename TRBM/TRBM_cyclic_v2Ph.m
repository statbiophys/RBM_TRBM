function [ Ph_l, mE_l ] = TRBM_cyclic_v2Ph( M, v )
%RBM_V2PH

% M : model
% v : vector column of visible units

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

Ph_l = TRBM_cyclic_v2mEh( M, v);
if nargout>1
    mE_l = Ph_l; % mE of the +1 state
end
Ph_l = sigmoid_approx( Ph_l );

end

