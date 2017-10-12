function [ Ph_l, mE_l ] = RBM_v2Ph( M, v )
%RBM_V2PH
% M : model
% v : vector column of visible units
% Ph_l : vector column of mean hidden states

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

Ph_l =  bsxfun(@plus, M.w*v, M.b' );
if nargout>1
    mE_l = Ph_l; % mE of the +1 state
end
Ph_l = sigmoid_approx( Ph_l );

end

