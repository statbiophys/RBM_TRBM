function [ logZ ] = RBM_logZ_Annealed_Discriminance_Sampling( M, K, Nt )
%RBM_Z_ANNEALED_IMPORTANCE_SAMPLING
% computes the partition function Z for RBM models using Annealed Discriminance Sampling
% the partition function Z is the normalization factor for the probability of 
% responses in single time bins, which log energy is computed
% using RBM_mF

% algo from Estimating the Partition Function by Discriminance Sampling, Liu et al, 2015
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

    function [f,gr] = fun_opt(q_new, q_old, x) % x = 1/c = Z_k / Z_k+1 compared to article
        t_old = 1./(q_old +x);
        t_new = 1./(q_new +x);
        f = q_old*t_old' + q_new*t_new' - length(q_new);  % 1 because the w are L1 normalized
        gr= - q_old*(t_old.^2)' - q_new*(t_new.^2)' ;
    end

Ni = length(M.a);
Nj = length(M.b);
if  Nj>Ni
    Minv = struct('a',M.b,'b',M.a,'w',M.w');
    logZ = RBM_Z_Annealed_Discriminance_Sampling( Minv, K, Nt );
    
else
    %% solver
    opts = optimoptions('fsolve','Display','none','Jacobian','on');
    
    %% independent model
    pi_l = M.pi_l;
    magn_l = log(pi_l./(1-pi_l))';
    
    M_old = struct('a',magn_l,'b',M.b*0,'w',M.w*0,'pi_l',M.pi_l); % RBM with independent neurons
    logZ =  sum(log(1./(1-pi_l))) + Nj*log(2); % logZ for p0
    
    v_l = bsxfun(@plus, rand(Ni, Nt), -pi_l)<0; % % generate independent responses
    %%
    r0_search = 1;
    for k = 1:K
        
        M_new = struct('a',magn_l*((K-k)/K) + M.a*(k/K),'b',M.b*(k/K),'w',M.w*(k/K),'pi_l',M.pi_l); % model at step k
        
        q_old = exp(RBM_mF(M_old,v_l) - RBM_mF(M_new,v_l)); % ratio of likelihoods for responses at step k-1
        v_l = RBM_simulate(M_new,v_l,1); % sampling of responses at step k
        q_new = exp(RBM_mF(M_old,v_l) - RBM_mF(M_new,v_l)); %  ratio of likelihoods for responses at step k
        
        r0_search =  fsolve(@(t) fun_opt(q_new, q_old, t), r0_search, opts);
        
        logZ = logZ - log(r0_search(end));
        step_notification(k,100,1000);
        
        M_old = M_new;
    end
    
end
end

