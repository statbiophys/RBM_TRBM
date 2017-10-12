function [ mF_l ] = TRBM_mF_cyclic( M, Nb_perseq, v_l, is_hidden)

% INPUT: 
% M: TRBM mode
% Nb_perseq: number of time bin per sequence
% v_l: responses: binary vector of size (Ni, Nseq*Nb_perseq)
% is_hidden: not used (exists only for homogeneity with RBM_mF_cyclic)

% OUTPUT: 
% mF_l: list of minus energy for each sequence, of size (1, Nseq)
% because the likelihood of the response in Nb_perseq time bins depends on responses at the boundary,
% which are not available, boundary responses are replaced by the firing rate of each neuron in the energy

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

[~, Nt] = size(v_l);
Nseq = Nt/Nb_perseq;
Nj = length(M.b);

Tmem = size(M.w,3);
v_d = TRBM_add_mean_field( v_l, Nb_perseq, Tmem, M.pi_l );

if nargin>=4 && is_hidden
    error('Case not available');
else
    switch 4 % see RBM_mF for details
        %         case 0
        %             xj_l = bsxfun(@plus, M.w*v_l, M.b');
        %             pj_l = sigmoid(xj_l);
        %             mF_l = M.a*v_l + sum( pj_l.*xj_l ,1) -sum( pj_l.*log(pj_l) + (1-pj_l).*log(1-pj_l) ,1);
        case 1
            [~, xj_l ] = TRBM_cyclic_v2Ph( M, v_d );
            mF_l1 = M.a*v_l;
            mF_l1 = sum(reshape(mF_l1,[Nb_perseq, Nseq]),1);
            xj_l = reshape(xj_l, [Nj*(Nb_perseq+Tmem-1), Nseq]);
            mF_l = mF_l1 + sum( log(1+exp(xj_l)) ,1);
        case 4
            mF_l1 = M.a*v_l;
            mF_l1 = sum(reshape(mF_l1,[Nb_perseq, Nseq]),1);
            
            xj_l = TRBM_cyclic_v2mEh( M, v_d );
            xj_l = reshape(xj_l, [Nj*(Nb_perseq+Tmem-1), Nseq]);
            
            mF_l = mF_l1 + sub_mF_prodlogexp_approx(xj_l, 5);
    end
end

