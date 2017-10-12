%% Load data
% two example responses
load(['~/Desktop/temp/response_data']);
resp1=resp_test(:,1:Nb_perseq);
resp2=resp_test(:,Nb_perseq+(1:Nb_perseq));

%  resp1 and resp2 are two responses of size (Nne, Nt) between which we
% compute the distance for the RBM and TRBM metrics. The number of time 
% bin Nt is not necessarily equal to Nb_perseq

% load RBM Model
M_RBM = load('~/Desktop/temp/response_data_RBM_Nj4_batch30_alpha0p01_beta1em05_moment1_b0m0p4_Nstep50_Nne60','M');
M_RBM = M_RBM.M;

% load TRBM Model
M_TRBM = load('~/Desktop/temp/response_data_TRBM_Tmem2_Nj4_batch10_alpha0p004_beta1em05_moment1_b0m0p4_Nstep100_Nne60');
M_TRBM = M_TRBM.M;
%% Euclidean RBM metric
Nne = size(resp1,1);

mh1 = RBM_v2Ph(M_RBM, resp1); % mean response of hidden units for resp1
mh2 = RBM_v2Ph(M_RBM, resp2);

d_RBM_Euclidean = func_dist_Ln( mh1, mh2, 2); % Eucliean RBM metric 

%% Semantic RBM metric
% one needs to estimate the covariance of responses predicted by the model
% in order to compute the distance kernel RBM_ker
% This involves generating multiple responses for the RBM. This takes time,
% but is only done once in order to define the kernel of the metric, and 
% is the same for all pairs of responses for which one wants to compute the
% distance
% Furthermode, the same kernel can be used to compute the distance between 
% pairs of responses of different lengths ( although responses in the same 
% pair should have same length )

T = 50000; % number of time bins used for simulation
resp_indep = bsxfun(@plus, rand(Nne, T), -M_RBM.pi_l)<0; % responses for independent model
resp_M = RBM_simulate( M_RBM, resp_indep, 300); % iterations of Block Gibbs sampling
RBM_ker = M_RBM.w*cov(resp_M')*(M_RBM.w'); % metric kernel

d_RBM_Semantic = dist_RBM_ker( mh1, mh2, RBM_ker ); % Semantic RBM metric

%% Euclidean TRBM metric
Tmem = size(M_TRBM.w, 3);

 X = repmat(M_TRBM.pi_l,[1 (Tmem-1)]); 
 mh1  = TRBM_cyclic_v2Ph(M_TRBM, [X, resp1]); 
 mh2  = TRBM_cyclic_v2Ph(M_TRBM, [X, resp2]); % mean hidden units, with mean field
 % approximation for responses out of boundaries

d_TRBM_Euclidean =  func_dist_Ln( mh1(:), mh2(:), 2);

%% Semantic TRBM metric
% similar to RBM but with temporal interactions.
% BEWARE: here the kernel TRBM_ker is specific to the length of responses 
% between which the distance is computed, in order to account for temporal 
% correlations

T = 50000; % number of time bins used for simulation
resp_indep = bsxfun(@plus, rand(Nne, T), -M_TRBM.pi_l)<0; % responses for independent model
resp_M = TRBM_simulate_cyclic( M_TRBM, resp_indep, 300); % iterations of Block Gibbs sampling
mE_l =  TRBM_cyclic_v2mEh( M_TRBM, resp_M);
TRBM_ker = cov(TRBM_unfold_time(mE_l, Inf, size(mh1,2))');
                    
d_TRBM_Semantic =  dist_RBM_ker( mh1(:), mh2(:), TRBM_ker);
