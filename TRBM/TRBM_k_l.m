function [ k_l, ck_l ] = TRBM_k_l( resp, Nb_perseq, Tcorr )
%  TRBM_k_l:
% k_l: list of the population count: number of spikes in the population in Tcorr consecutive time bins
% ck_l: distribution of the population count, with ck_l[i] the frequency of value (i-1) 

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

k_l = full(sum(resp,1));
k_l =  sum(TRBM_unfold_time(k_l, Nb_perseq, Tcorr),1);

Ni = size(resp,1);
ck_l = histc(k_l,0:(Ni*Tcorr));
ck_l = ck_l/sum(ck_l);

end

