function [ mF_l ] = RBM_mF( M, v_l, is_hidden)
% list of minus free energy
% if is_hidden = 1, input v is actually h

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

if nargin<3
    is_hidden=0;
end
if is_hidden
    Minv = struct('a',M.b,'b',M.a,'w',M.w');
    [ mF_l ] = RBM_mF( Minv, v_l, 0);
else
    switch 4
        case 0
            xj_l = bsxfun(@plus, M.w*v_l, M.b');
            pj_l = sigmoid(xj_l);
            mF_l = M.a*v_l + sum( pj_l.*xj_l ,1) -sum( pj_l.*log(pj_l) + (1-pj_l).*log(1-pj_l) ,1);
        case 1
            xj_l = bsxfun(@plus, M.w*v_l, M.b');
            mF_l = M.a*v_l + sum( log(1+exp(xj_l)) ,1);
        case 2 % fast but goes to +Inf easily
            xj_l = bsxfun(@plus, M.w*v_l, M.b');
            mF_l = M.a*v_l +  log(prod(1+exp(xj_l),1));
        case 3
            mF_l = sub_RBM_mF_mex(M,M.w*v_l) + M.a*v_l;
        case 4
            mF_l = sub_RBM_mF_mex_approx(M,M.w*v_l, 5) + M.a*v_l;
        case 5 
            xj_l = bsxfun(@plus, M.w*v_l, M.b');
            x_max = 5;
            ke = xj_l < x_max;
 
            exj_l = ones(size(xj_l));
            exj_l(ke) = 1+exp(xj_l(ke));
            
            xj_l(ke) = 0;
             
            mF_l = M.a*v_l + sum(xj_l, 1) + log(prod(exj_l,1));
            if any(isinf(mF_l(:)))
                error('Inf found');
            end
    end

    
% if 0
%     %% test
%     Ni = 50;
%     Nj = 100;
%     Nt = 10000;
%     
%     M.a = ones(1,Ni); %normrnd(0,1,1,Ni)*3;
%     M.b = normrnd(0,1,1,Nj)*3;
%     M.w = normrnd(0,1,Nj,Ni)*3;
%     
%     v_l = rand(Ni, Nt)>0.5;
%     %%
%     RBM_mF( M, v_l )
%     %%
%     tic
%     for u = 1:100
%         [ mF_l ] = RBM_mF( M, v_l, 0);
%     end
%     toc
% end
end

