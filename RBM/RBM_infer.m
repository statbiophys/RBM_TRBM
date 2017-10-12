function [ M, testvals, L ] = RBM_infer( v_l, batch_size, Nbin_CD, Nstep, Nstep_CD, param, v_l_heldout, M0)
%RBM_TRAIN
% INPUT:
% v_l : response: vector of size (Nneuron, Nt)
% batch_size: size of batches used (in number of time bins)
% Nb_CD: number of sequences used for Contrastive Divergence
% Nstep: number of steps during inference
% Nstep_CD: number of Block Gibbs Sampling steps during Contrastive Divergence
% param: structure containing parameters
% v_l_heldout: responses in testing set. Used only to compute energy or likelihood of 
%              responses for intermediate models during learning. Can be "[]" is not used
% M0 : initial model (optional)

% OUPUT:
% M: model
% testvals: structure containing some statistics computed for intermediate models (specified in 'param')
% L: some other statistics computed for intermediate models (modify functin code below, convenient for debugging)

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

    function ke_il = fun_ke_il( N, k)
        if k == N
            ke_il = 1:N;
        else
            %         if 0
            %             ke_il = randperm(N);
            %             ke_il = sort(ke_il(1:k),'ascend');
            %         else % faster
            ke_il = zeros(1,N);
            ke_il(1:k) = 1;
            ke_il = find(ke_il(randperm(N)));
            %         end
        end
    end

%% check inputs
if full(any( (v_l(:) ~=0) & ( (v_l(:) ~=1) ))) ||  full(any( (v_l_heldout(:) ~=0) & ( (v_l_heldout(:) ~=1) )))
    error('v_l and v_l_heldout must be binary');
end

if Nstep_CD.comp(1) ~= 0
    error('Nstep_CD.comp(1) must be 0');
end

Nj = param.Nj;
fprintf(['Learning RBM with Nj: ' int2str(Nj) ' hidden units per time bin\n']);

[Ni, Nt] = size(v_l);
Nbatch = ceil(Nt/batch_size);

%% initialization : see Hinton 2010 "A practical guide to training RBM"
if exist('M0','var')
    M = M0;
else
    M = struct();
end
M.pi_l = full(mean(v_l,2)); % Mean response per neuron. Convenient for mean-field at boundaries
if ~isfield(M,'a');
    M.a = log(M.pi_l./(1-M.pi_l))';
end
if ~isfield(M,'b');
    M.b = zeros(1,Nj);
end
if ~isfield(M,'w');
    M.w = normrnd(0,0.01,Nj,Ni);
end

%% v_CD: responses used for Contrastive Divergence
kei = fun_ke_il( Nt, Nbin_CD);
v_CD = v_l(:,kei);

%% Parameters
if param.is_moment
    Da_old = 0;
    Db_old = 0;
    Dw_old = 0;
end

alpha_a = param.alpha;
alpha_b = param.alpha;
alpha_w = param.alpha;

if nargout>1 % used to output testvals and L, showing some statistics for models during inference
    testvals = struct();
    Ncomp = Nstep*Nbatch;
    if param.Ncomp_mF % negative free Energy mF (m stands for 'minus')
        testvals.mF_l = zeros(2,floor(Ncomp/param.Ncomp_mF));
    end
    if param.Ncomp_Z % partition function Z
        testvals.Z_l = zeros(1,floor(Ncomp/param.Ncomp_Z));
    end
    if param.Ncomp_ph % mean response of hidden units
        testvals.ph = zeros(Nj,floor(Ncomp/param.Ncomp_ph));
    end
    if param.Ncomp_ri % mean response of neurons
        testvals.ri = zeros(Ni,floor(Ncomp/param.Ncomp_ri));
    end
    if param.Ncomp_Cij % pearson correlation between pairs of neurons
        testvals.Cij = zeros(Ni*(Ni-1)/2,floor(Ncomp/param.Ncomp_Cij));
    end    
    if nargout>2
        L = [];
    end
end

%% Mean observables before training
if param.Ncomp_ph
    [ ~,~, pj_l_mod ] = RBM_simulate(M, v_CD(:,mod(1:param.ph_Nt, end)+1), param.ph_Nstep );
    testvals.ph(:,1) =  pj_l_mod;
end

%% Training
fprintf(['Sarting for ' int2str(Nstep) ' steps with ' int2str(Nbatch) ' batches\n']);

comp_i = 0;
for step_i = 1:Nstep
    ord = randperm(Nt);
    for batch_i = 1:Nbatch
        comp_i = comp_i+1;
        
        v_p = v_l(:,ord(((batch_i-1)*batch_size+1):min((batch_i*batch_size),end)));
        v_p = full(v_p); % sub-part of responses selected for computations ( '_p' stands for 'partial')
                
        Nstep_CD_t = Nstep_CD.val(find(comp_i>=Nstep_CD.comp,1,'last'));
        [ v_CD, pi_l_mod, pj_l_mod, pji_l_mod ] = RBM_simulate( M, v_CD, Nstep_CD_t );
                
        da = mean(v_p,2)' - pi_l_mod';
        
        pj_l = RBM_v2Ph(M, v_p);
        
        pji_l = pj_l*v_p';
        pji_l = pji_l/size(v_p,2);        
        dw = pji_l - pji_l_mod;
        
        %         if 0
        %             pj_l = mean(pj_l,2);
        %             db = pj_l' - pj_l_mod';
        %         else % suggested by Hinton
        db =  mean(pj_l>rand(size(pj_l)),2)' - pj_l_mod';
        %         end
        
        if nargout>2  % un-comment what one wants to compute
            l = [];
            %                     l = [l; norm(M.a,1); norm(M.b,1); norm(M.w(:),1)];
            %                     l = [l; norm(da,1); norm(db,1); norm(dw(:),1)];
            %                     l = [l; norm(da,1)/norm(M.a,1); norm(db,1)/norm(M.b,1); norm(dw(:),1)/norm(M.w(:),1)];
            %                     l = [l; M.a(:); M.b(:); M.w(:)];
            %                     l = [l; alpha_a; alpha_b; alpha_w];
            L = [L, l];
        end
        
        %% update model
        if ~param.is_moment
            M.a = M.a + alpha_a*da;
            M.b = M.b + alpha_b*db;
            M.w = M.w*(1-alpha_w*param.beta) + alpha_w*dw;
        else
            mom = param.moment_f(step_i);
            Da = alpha_a*da + mom*Da_old;
            Db = alpha_b*db + mom*Db_old;
            Dw = alpha_w*(-param.beta*M.w + dw) + mom*Dw_old;
            
            M.a = M.a + Da;
            M.b = M.b + Db;
            M.w = M.w + Dw;
            
            Da_old = Da;
            Db_old = Db;
            Dw_old = Dw;
        end
        
        %%
        if nargout>1 % compute some statistics (optional)
            if ~mod(comp_i,param.Ncomp_mF)
                testvals.mF_l(:,comp_i/param.Ncomp_mF) = [mean(RBM_mF(M,v_l)); mean(RBM_mF(M,v_l_heldout))];
            end
            if ~mod(comp_i,param.Ncomp_Z)
                if Nj<25
                    testvals.Z_l(1,comp_i/param.Ncomp_Z) = RBM_Z_exact(M);
                else
                    [testvals.Z_l(1,comp_i/param.Ncomp_Z)] = ... 
                       exp(RBM_logZ_Annealed_Importance_Sampling(M, param.AIS_Nstep, param.AIS_Nt));
                    %                     [testvals.Z_l(1,comp_i/param.Ncomp_Z)] = ...
                    %                         exp(RBM_logZ_Annealed_Discriminance_Sampling(M, param.ADS_Nstep, param.ADS_Nt));
                end
            end           
            if ~mod(comp_i,param.Ncomp_ph)
                [ ~,~, pj_l_mod ] = RBM_simulate(M, v_CD(:,mod(1:param.ph_Nt,end)+1), param.ph_Nstep );
                testvals.ph(:,comp_i/param.Ncomp_ph +1) =  pj_l_mod;
            end
            if ~mod(comp_i,param.Ncomp_ri)
                testvals.ri(:,comp_i/param.Ncomp_ri) = pi_l_mod;
            end
            if ~mod(comp_i,param.Ncomp_Cij)
                testvals.Cij(:,comp_i/param.Ncomp_Cij) = triu_l(corr(v_CD'));
            end
        end
%         step_notification(batch_i, 10, 100); % un-comment to display batch numbers
    end
    fprintf(['done step ' int2str(step_i) '\n']);
end

if nargout>1 && param.Ncomp_Z && param.Ncomp_mF && ~isempty(param.Ncomp_Z) && isequal(size(testvals.Z_l,2), size(testvals.mF_l,2))
    testvals.loglike = bsxfun(@plus, testvals.mF_l, - log(testvals.Z_l));
end

fprintf('\n RBM inference ended\n');
end

