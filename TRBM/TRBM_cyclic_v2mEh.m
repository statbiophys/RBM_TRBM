function [ mE_l ] = TRBM_cyclic_v2mEh( M, v )
% M : model
% v : vector column of visible units

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

[Nj, ~, Tmem] = size(M.w);

tm = nmodeproduct(M.w, v', 2, 1);
kei = (0:(Nj-1))*Tmem;
t2 = tm(:,kei+1);
if ~verLessThan('matlab', '8.3')
    for d = 1:(Tmem-1)
        t2 = t2 + circshift(tm(:, d +1 + kei),d, 1);
    end
else
    for d = 1:(Tmem-1)
        t2 = t2 + circshift_asOLD(tm(:, d +1 + kei),d, 1);
    end
end

mE_l =  bsxfun(@plus, t2', M.b' );

end

