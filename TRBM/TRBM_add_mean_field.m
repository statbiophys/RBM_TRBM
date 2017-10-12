function [ v2 ] = TRBM_add_mean_field( v_l, Nb_perseq, Tmem, pi_l )
%TRBM_ADD_MEAN_FIELD 

% v_l is a series of [ Ni,  Nseq*Nb_perseq]
% this function inserts (Tmem-1) mean field approximations of v (i.e. its mean, pi_l) between
% each sequence

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

[Ni, Nt] = size(v_l);
Nseq = Nt/Nb_perseq;

v2 = reshape( v_l, [Ni*Nb_perseq, Nseq]);

v2 = [v2; repmat( pi_l, [(Tmem-1), Nseq] )];
v2 = reshape( v2, [Ni, Nseq*(Nb_perseq + Tmem -1)]);

if 0
    %%      test
    Nb_perseq = 4;
    Tmem = 4;
    Nseq = 2;
    Ni = 3;
    pi_l = 100 + (1:Ni)'
    v_l = rand(Ni, Nseq*Nb_perseq)
    
    TRBM_add_mean_field( v_l, Nb_perseq, Tmem, pi_l )
end

