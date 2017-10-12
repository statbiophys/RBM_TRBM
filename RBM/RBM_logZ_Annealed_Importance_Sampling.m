function [ logZ ] = RBM_logZ_Annealed_Importance_Sampling( M, K, Nt)
%RBM_Z_ANNEALED_IMPORTANCE_SAMPLING
% computes the partition function Z for RBM models using Annealed Importance Sampling
% the partition function Z is the normalization factor for the probability of 
% responses in single time bins, which log energy is computed
% using RBM_mF
%
% algo from Learning and evaluating Boltzmann machines, Salakhutdinov, 2008
% with E_k = (k/K) * E_M + (K-k)/k * E0
% with E0 corresponding to a model of independent sigma with approx same firing rate

% INPUT:
% M: RBM model
% M must contain M.pi_l, column vector of mean response for each neuron
% K: number of temperatures used
% Nt: number of time bins generated for each temperature

% OUTPUT
% logZ: logarithm of the partition function Z

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

[Nj, Ni] = size(M.w);

if  Nj>Ni
    Minv = struct('a',M.b,'b',M.a,'w',M.w');
    [ logZ ] = RBM_logZperbin_Annealed_Importance_Sampling( Minv, K, Nt);
    
else
    %%
    %     pi_l = mean( RBM_simulate(M, rand(Ni, Nt)>0.5, 50), 2); 
    %     pi_l(pi_l == 0) = min(pi_l(pi_l>0)); % avoid degenerate model
    pi_l = M.pi_l;
    
    magn_l = log(pi_l./(1-pi_l))';
    
    %%
    v_l = bsxfun(@plus, rand(Ni, Nt), -pi_l)<0; % p0 : independent model
    
    M_t = struct('a',magn_l,'b',M.b*0,'w',M.w*0,'pi_l',M.pi_l);
       
    logZ =  sum(log(1./(1-pi_l))) + Nj*log(2); % logZ for p0
    
    for k = 1:K
        mF_old = RBM_mF( M_t, v_l); % old model
        M_t = struct('a',magn_l*((K-k)/K) + M.a*(k/K),'b',M.b*(k/K),'w',M.w*(k/K),'pi_l',M.pi_l); % model at step k
        mF_new = RBM_mF( M_t, v_l); % new model
        logZ = logZ + log(mean(exp(mF_new - mF_old)));
        v_l = RBM_simulate(M_t, v_l, 1); % Block Gibbs Sampling
        
        step_notification(k, 100, 1000);
    end
end
end

