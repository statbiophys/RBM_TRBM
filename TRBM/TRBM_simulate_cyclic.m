function [ v_l, pi_l, pj_l, pji_m ] = TRBM_simulate_cyclic( M, v_l, N )
% TRBM_SIMULATE
% simulates TRBM during N time steps
% using cyclic boundary conditions

% INPUT
% M: model
% v_l: vector or responses for initialization
% number of iterations

% OUTPUT
% v_l: vector of responses
% pi_l: list of firing rates of visible units
% pj_l: list of firing rates of hidden units
% pji_m: tensor of linear correlation between visible units and hidden units: m(j,i,d) = < h_(j,t+d-1) sigma_(i,t) > 

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

[Ni, Nt] = size(v_l);
Nj = size(M.w,1);

for n = 1:(N-1)
    h_l = TRBM_cyclic_v2Ph(M, v_l)>rand(Nj, Nt);
    v_l = TRBM_cyclic_h2Pv(M, h_l)>rand(Ni, Nt);
    %     step_notification(n, 1, 10);
end

h_l = TRBM_cyclic_v2Ph(M, v_l)>rand(Nj, Nt);

if nargout>1
    pi_l = TRBM_cyclic_h2Pv(M, h_l);
    v_l = pi_l>rand(Ni, Nt);
    
    pj_l = TRBM_cyclic_v2Ph(M, v_l);
    
    pji_m = zeros(size(M.w));
    if ~verLessThan('matlab', '8.3')
        for d = 1:size(M.w,3)
            pji_m(:,:,d) = circshift(pj_l,1-d,2)*v_l';
        end
    else
        for d = 1:size(M.w,3)
            pji_m(:,:,d) = circshift_asOLD(pj_l,1-d,2)*v_l';
        end
    end
    
    pji_m = pji_m/Nt;
    
    pi_l = mean(pi_l,2);
    pj_l = mean(pj_l,2);
else
    v_l = TRBM_cyclic_h2Pv(M, h_l)>rand(Ni, Nt);
end
end

