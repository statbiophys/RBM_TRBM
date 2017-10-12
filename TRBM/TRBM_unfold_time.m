function V = TRBM_unfold_time( v, Nb_perseq, B)
%TRBM_UNFOLD_TIME
% time spatialization for B consective time bins. 

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

if isinf(Nb_perseq)
    Nb_perseq = size(v,2);
end

V = v;
if ~verLessThan('matlab', '8.3')
    for b = 1:(B-1)
        V = [V; circshift(v, -b,2)];
    end
else
    for b = 1:(B-1)
        V = [V; circshift_asOLD(v, -b,2)];
    end
end
ke = 1:size(v,2);
ke = mod(ke-1,Nb_perseq) <= (Nb_perseq - B +1) -1;
V = V(:, ke);

% if 0
%     %% test
%     Nb_perseq = 15;
%     Nseq = 8;
%     Ni = 2;
%     v = rand(Ni, Nb_perseq*Nseq);
%
%     B = 3;
%     %%
%     tic
%     V = TRBM_unfold_time( v, Nb_perseq, B);
%     toc
%     %%
%     tic
%     V2 =  TRBM_unfold_time( v, Nb_perseq, B);
%     toc
%     %%
%     isequal(V,V2)
%     %%
%     v(:,1:10)
%     V(:,1:10)
%     V2(:,1:10)
% end
end