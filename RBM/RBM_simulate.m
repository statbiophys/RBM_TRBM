function [ v_l, pi_l, pj_l, pji_l ] = RBM_simulate( M, v_l, N )
%RBM_SIMULATE
% simulates RBM during N time steps

% INPUT
% M: model
% v_l: vector or responses for initialization

% OUTPUT
% v_l: vector of responses
% pi_l: list of firing rates of visible units
% pj_l: list of firing rates of hidden units
% pji_l: matrix of correlation between visible units and hidden units: m(j,i) = < h_j sigma_i > 

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

Nt = size(v_l,2);
Nj = size(M.w,1);

for n = 1:(N-1)
    h_l = RBM_v2Ph(M, v_l)>rand(Nj, Nt);
    v_l = RBM_h2Pv(M, h_l)>rand(size(v_l));
end

h_l = RBM_v2Ph(M, v_l)>rand(Nj, Nt);

if nargout>1
    pi_l = RBM_h2Pv(M, h_l);
    v_l = pi_l>rand(size(v_l));

    pj_l = RBM_v2Ph(M, v_l);
    
    pji_l = pj_l*v_l';
    pji_l=pji_l/Nt;
    
    pi_l = mean(pi_l,2);
    pj_l = mean(pj_l,2);
else
    v_l = RBM_h2Pv(M, h_l)>rand(size(v_l)); 
end
end

