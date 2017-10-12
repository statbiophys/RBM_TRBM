function [ logZ ] = TRBM_logZ_Annealed_Importance_Sampling( M, K, Nt, Nb_perseq)
%RBM_Z_ANNEALED_IMPORTANCE_SAMPLING
% computes the partition function Z for TRBM models using Annealed Importance Sampling
% the partition function Z is the normalization factor for the probability of 
% responses in Nb_perseq consecutive time bins, which log energy is computed
% using TRBM_mF_cyclic
% BEWARE: here Z is specific to the length of responses Nb_perseq !
%
% algo from Learning and evaluating Boltzmann machines, Salakhutdinov, 2008
% with E_k = (k/K) * E_M + (K-k)/k * E0
% with E0 corresponding to a model of independent sigma with approx same firing rate

% INPUT:
% M: TRBM model
% M must contain M.pi_l, column vector of mean response for each neuron
% K: number of temperatures used
% Nt: number of time bins generated for each temperature
% Nb_perseq: length (in time bins) for responses considered

% OUTPUT
% logZ: logarithm of the partition function Z

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

Nt = Nb_perseq*ceil(Nt/Nb_perseq); % make sure Nb_perseq divides Nt

[Nj, Ni, Tmem] = size(M.w);

if  Nj>Ni
    Minv = struct('a',M.b,'b',M.a,'w',permute(M.w,[2 1 3]));
    [ Z ] = TRBM_logZperbin_Annealed_Importance_Sampling( Minv, K, Nt);
    
else
    %%
    pi_l = M.pi_l;
    if any(pi_l==0); error('Elements of M.pi_l must be positive'); end;
    magn_l = log(pi_l./(1-pi_l))';
    
    %%
    v_l = bsxfun(@plus, rand(Ni, Nt), -pi_l)<0; % generate independent responses
    
    M_t = struct('a',magn_l,'b',M.b*0,'w',M.w*0,'pi_l',M.pi_l); % p0 : independent model 
    logZ =  Nb_perseq*sum(log(1./(1-pi_l))) + Nj*(Nb_perseq+Tmem-1)*log(2); % logZ for p0
    
    for k = 1:K
        mF_old = TRBM_mF_cyclic( M_t, Nb_perseq, v_l); % old model
        M_t = struct('a',magn_l*((K-k)/K) + M.a*(k/K),'b',M.b*(k/K),'w',M.w*(k/K),'pi_l',M.pi_l); % model at step k
        mF_new = TRBM_mF_cyclic( M_t, Nb_perseq, v_l); % new model
        logZ = logZ + log(mean(exp(mF_new - mF_old)));
        v_l = TRBM_simulate_cyclic(M_t, v_l, 1); % Block Gibbs Sampling
        
        step_notification(k, 100, 1000);
    end
end
end

