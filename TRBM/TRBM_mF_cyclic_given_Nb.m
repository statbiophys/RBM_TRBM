function [ mF_l ] = TRBM_mF_cyclic_given_Nb( M, Nb_perseq, v_l, loglikel_Nbin)
% Computes mean free energy for responses of size loglikel_Nbin consecutive time bins

% INPUT: 
% M: TRBM mode
% Nb_perseq: number of time bin per sequence
% v_l: responses: binary vector of size (Ni, Nseq*Nb_perseq)

% OUTPUT: 
% mF_l: list of minus energy for each sequence, of size (1, Nseq)
% because the likelihood of the response in Nb_perseq time bins depends on responses at the boundary,
% which are not available, boundary responses are replaced by the firing rate of each neuron in the energy

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

% Creates responses of loglikel_Nbin consecutive time bins
ke = zeros(1,Nb_perseq);
ke(1:loglikel_Nbin*floor(Nb_perseq/loglikel_Nbin))=1;
v_l = v_l(:,repmat(ke,[1 size(v_l,2)/Nb_perseq])>0);

% Computes mF for responses in loglikel_Nbin consective time bins
[ mF_l ] = TRBM_mF_cyclic( M, loglikel_Nbin, v_l);
end

