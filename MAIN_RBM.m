%% Inference
%% Data
fpath_data = ['~/Desktop/temp/response_data']; % path to data. Must contain variables:
% 'resp_train' : binary sparse matrix of size (Nneuron, Nbin) with  0 (silence) and 1 (spikes)

fpath_result_RBM = ['~/Desktop/temp/response_data_RBM']; % path of file where results are saved

load(fpath_data);
if 0 % restrict to some subpopulation only (optionnal)
    kei_i = 1:10 % list of indices of chosen neurons
    resp_train = resp_train(kei_i,:);
    resp_test = resp_test(kei_i,:);
end

resp_train = resp_train>0; % make sure resp are binary
resp_test = resp_test>0;

[Nne,Nt_train] = size(resp_train);

%% Parameters
param.Nj = 20;  % number of hidden units

Nstep = 100; % number of iterations of the inference
batch_size = 10; % size (in number of sequences) of the batches used for inference. Each iteration 
% of the inference treates Nbatch ~= Nt/batch_size  batches one after the other.
% Treating one batch consists of one 'computation', and is indexed by comp_i = 1,...,Nbatch*Nstep
% in the algorithm.

Nbatch = ceil(Nt_train/batch_size);

Nbin_CD = 500; % number of responses generated during contrastive divergence (CD)

Nstep_CD.val = 1:4; % number of Block Gibbs sampling steps during CD.
Nstep_CD.comp = (0:3)*40*Nbatch; % At computation comp_i, the algorithm finds the largest k such 
% that Nstep_CD.comp(k) <= comp_i, and then CD is done using Nstep_CD.val(k) steps.
% Example: Nstep_CD.val = 1:4 ,  Nstep_CD.comp = (0:3)*40*Nbatch;   % the number of iterations
% increases from 1 to 4
% Other example:  Nstep_CD.val = 4; , Nstep_CD.comp = 0;   % to always have 4 iterations

param.alpha = 0.01; % learning coefficient
param.beta = 10.^-5; % L2 regularization

b0 = -0.4; % inital value for fields b_j

param.is_moment = 1; % if 1, uses moments from the gradient. Then, at computation
% comp_i, the scaling coefficient of the moment is param.moment_f(comp_i)
comp_switch = Nbatch*10;
param.moment_f = @(x) (0.5*(x<comp_switch) + 0.9*(x>=comp_switch));

%% additional parameters
% Some statistics can be computed for models during learning, and are stored in the structure 'testvals'
% In order to save time, statistics are not computed at each computations, but only every 'Ncomp' computations.
% This Ncomp is different for each statistic, and are set here.
% If Ncomp is set to 0, the corresponding statistic is not computed.
% These can take a significant amount of time: skip it once the learning parameters are set.

param.Ncomp_ri = 0*Nbatch; % firing rate of neurons
param.Ncomp_Cij = 0*Nbatch; % Pearson correlation between all pairs of neurons

param.Ncomp_mF = 0*Nbatch; % negative mean energy of responses in training and testing set (separately)

param.Ncomp_Z = 0*Nbatch; % partition function
param.AIS_Nstep = 5000;  % number of steps used for computation of Z with Annealed Importance sampling
param.AIS_Nt = 5000;  % number of time bins used for computation of Z with Annealed Importance sampling

param.Ncomp_ph = 0*Nbatch; %Nbatch; % firing rate of hidden units. Computed with a simulation with parameters:
param.ph_Nt = 10000; % number of time bins used for computation of ph
param.ph_Nstep = 15; % number of Block Gibbs Sampling iterations used for computation of ph

%% Initialization and inference
M0 = struct('b',b0*ones(1,param.Nj));

tic
[ M, testvals] = RBM_infer( resp_train, batch_size, Nbin_CD,  Nstep, Nstep_CD, param, resp_test, M0);
toc
comp_time = toc;

fprintf(['\n -----  done RBM Nj: ' int2str(param.Nj)  ', beta: ' num2str(param.beta) '\n']);

if 1
    %%
    fpath_result_RBM_wparams = [fpath_result_RBM '_Nj' int2str(param.Nj) ...
        '_batch' int2str(batch_size) '_alpha' num2str_dot2p(param.alpha) ...
        '_beta' num2str_dot2p(param.beta) '_moment'  int2str(param.is_moment) ...
        '_b0' num2str_dot2p(b0) ...
        '_Nstep' int2str(Nstep) '_Nne' int2str(Nne)]; % path to save results, with parameters
    
    save(fpath_result_RBM_wparams,'M','testvals','param','comp_time');
end

if 0
    %% Plots: learning
    % these plots are available only if these values have been computed 
    % during the exeriment. This is achieved by setting param.Ncomp_ ... parameters
    % with positive value
    %%
    f=figure;
    plot( testvals.mF_l' );
    
    m_l= mean(testvals.ph,1);
    std_l = std(testvals.ph);
    
    legend({'train','test'},'Location','NorthWest');
    xlabel('step');
    ylabel('free Energy');
    %%
    f=figure;
    plot( testvals.loglike' );

    legend({'train','test'},'Location','NorthWest');
    xlabel('step');
    ylabel('mean loglikelihood');
    %%
    f=figure; hold on
    plot(m_l,'color','k');
    plot(m_l - std_l,'color','r');
    plot(m_l + std_l,'color','r');
    
    xlabel('step');
    ylabel('P(h_j = 1): mean ± std');
    %%
    f=figure; hold on
    ri_err = bsxfun(@plus,testvals.ri,-M.pi_l);
    ri_err_m = mean(ri_err,1)';
    ri_err_std = std(ri_err,0,1)';
    plot(ri_err');
    plot(ri_err_m,'r','LineWidth',2);
    plot(ri_err_m+ri_err_std,'r','LineWidth',2);
    plot(ri_err_m-ri_err_std,'r','LineWidth',2);
    plot_hline(0);
    
    xlabel('step');
    ylabel('error in P(h_j = 1): mean ± std');
    %%
    f=figure; hold on
    Cij_err = bsxfun(@plus,testvals.Cij,-triu_l(corr(resp_test')));
    Cij_err_m = mean(Cij_err,1)';
    Cij_err_std = std(Cij_err,0,1)';
    plot(Cij_err');
    plot(Cij_err_m,'r','LineWidth',2);
    plot(Cij_err_m+Cij_err_std,'r','LineWidth',2);
    plot(Cij_err_m-Cij_err_std,'r','LineWidth',2);
    plot_hline(0);
    
    xlabel('step');
    ylabel('error in corr(\sigma_i, \sigma_i): mean ± std');
end


%% Model quality
Tcorr_max=1; % code common to TRBM

r_test = mean(resp_test,2);
corr_test_m =  TRBM_corr_wdelay_l( resp_test, Nb_perseq, Tcorr_max );
ck_test_c = cell(1,Tcorr_max);
for t = 1:Tcorr_max
    [ ~, ck_test_c{t} ] = TRBM_k_l( resp_test, Nb_perseq, t);
end

r_train = mean(resp_train,2);
resp_indep = bsxfun(@plus, rand(size(resp_train)), -r_train)<0;
ck_indep_c = cell(1,Tcorr_max);
for t = 1:Tcorr_max
    [ ~, ck_indep_c{t} ] = TRBM_k_l( resp_indep, Inf, t);
end

resp_M =  RBM_simulate( M, resp_indep, 100);
r_l = mean(resp_M,2);
corr_m = corr_wdelay_l( resp_M, Tcorr_max );
ck_c = cell(1,Tcorr_max);
for t = 1:Tcorr_max
    [ ~, ck_c{t} ] = TRBM_k_l( resp_M, Inf, t );
end

if 1
    %% Plots: model quality
    %% Parameters for plots
    stim_rate = 50; % time bin frequency
    col_RBM =  [0.87 0.58 0.20];
    col_indep = [0.65 0.65 0.65];
    
    %% ****** r_i
    f=figure; hold on
    plot( r_test*stim_rate, r_l*stim_rate,'.','Color',col_RBM,'MarkerSize',12);
    
    resize_max_square;
    plot_identity_line;
    
    xlabel('data');
    ylabel('RBM');
    title('firing rates (Hz)');
    
    %% ****** P(K) (for figure: Tcorr = 1 and 5)
    Tcorr = 1;
    f=figure; hold on; ax = gca;
    
    Ymin = max([min(ck_test_c{Tcorr}), min(ck_c{Tcorr}), min(ck_indep_c{Tcorr})]);
    
    x_l = 0:Nne*Tcorr;
    ck_test_ste = sqrt(ck_test_c{Tcorr}.*(1-ck_test_c{Tcorr}))/sqrt(size(resp_test,2));
    plot_area_around_line(x_l, ck_test_c{Tcorr}, ck_test_ste,0.8*[1 1 1],1,Ymin);
    
    
    plot(x_l, ck_test_c{Tcorr},'k');
    plot(x_l, ck_indep_c{Tcorr},'Color',col_indep,'LineWidth',1.5);
    leg_c = {'','test data','independent'};
    
    plot(x_l, ck_c{Tcorr},':','LineWidth',3,'Color',col_RBM);
    leg_c{end+1} = ['RBM Nj: ' int2str(param.Nj)];
    
    ax.YScale='log';
    legend(leg_c,'Location','NorthEast'); legend boxoff;
    
    title(['Pop rate. Tcorr: ' int2str(Tcorr)]);
    xlabel('K');
    ylabel('probability');
    
    %% ******  C_ii', RBM
    f=figure; hold on
    ke = corr_test_m(:,1)<1;
    plot( corr_test_m(ke,1), corr_m(ke,1),'.','Color',col_RBM,'MarkerSize',12);
    
    resize_max_square;
    plot_identity_line;
    
    title(['Pairwise corr, corr: ' num2str(corr(corr_test_m(ke,1), corr_m(ke,1)),'%.3f')]);
    xlabel('data');
    ylabel('RBM');

end

%% Mean logLikelihood
% Compute partition function Z 
% [ logZ ] = log(RBM_Z_exact( M ));

Z_K = 5000; % number of temperatures
Z_Nt = 5000; % number of time bins for simuation
[ logZ ] = RBM_logZ_Annealed_Importance_Sampling( M, Z_K, Z_Nt); % log partition function with AIS algorithm
% [ logZ ] = RBM_logZ_Annealed_Discriminance_Sampling( M, Z_K, Z_Nt); % log partition function with ADS algorithm
% it seemed that AIS gives better estimates in some cases

% Mean response log likelihood
loglikel_train_m = mean(RBM_mF(M, resp_train)) - logZ
loglikel_test_m  = mean(RBM_mF(M, resp_test))  - logZ



