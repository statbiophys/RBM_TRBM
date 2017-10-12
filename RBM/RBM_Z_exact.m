function [ Z ] = RBM_Z_exact( M )
%RBM_LOGZ_EXACT
% Computes partition function exactly
%
% BEWARE: only possible for small number of visible units Ni or hidden units Nj 

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

Nj = length(M.b);
Ni = length(M.a);

if Nj>Ni
    Minv = struct('a',M.b,'b',M.a,'w',M.w');
    Z = RBM_Z_exact( Minv );
else
    if Nj>25
        error('Nj and Ni too large for exact computation');
    end
    v_l = permn([0 1],Nj)';
    [ mF_l ] = RBM_mF( M, v_l, 1);
    Z = sum(exp(mF_l));
end


% if 0
%     %% test
%     for i = 1:25
%         fprintf([int2str(i) ' : ']);
%         tic
%         v_l = permn([0 1],i);
%         toc
%     end
%     %%
%     Ni = 3;
%     Nj = 2;
%     Nt = 100000;
%     
%     M.a = ones(1,Ni); %normrnd(0,1,1,Ni)*3;
%     M.b = normrnd(0,1,1,Nj)*03;
%     M.w = normrnd(0,1,Nj,Ni)*0.3;
%     
%     v_l = rand(Ni, Nt)>0.5;
% 
%     %%
%     [ v_l] = RBM_simulate( M, v_l, 100 );
%  
%     [v_u,~,IC] = unique(v_l','rows'); v_u = v_u';
%     c = 0;
%     for i=1:size(v_u,2)
%         c(i) = mean(IC==i);
%     end
%     v_u
%     c
%     
%     c2 = exp(RBM_mF(M, v_u))/RBM_Z_exact( M )
% end
end

