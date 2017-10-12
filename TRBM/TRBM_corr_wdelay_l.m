function [ corr_l ] = TRBM_corr_wdelay_l( x, Nb_perseq, Tcorr )
%TRBM_CORR_WDELAY_L
% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

pi_l = mean(x,2);

%[Ni, Nt] = size(x);

% x = bsxfun(@plus, x, -pi_l(:));

x = TRBM_add_mean_field(  x, Nb_perseq, Tcorr, pi_l );

corr_l = corr_wdelay_l( x, Tcorr );

for t = 0:(Tcorr-1)
    corr_l(:,t+1) =  corr_l(:,t+1)*((Nb_perseq)/(Nb_perseq - t));
    step_notification(t+1, 1, 15);
end

% if 0
%     %% test
%     Tcorr = 3
%     Ni = 3;
%     Nt = 10000000;
%     x = normrnd(1,1,Ni,Nt);
%     x = [1 2 3; 0 1 -2; 0 0 1]*x;
%     x = circshift(x,1,2) - x;
%     x = circshift(x,1,2) + x;
%     C0 =  corr_wdelay_l(x, Tcorr)
%     
%     Nb_perseq = 10;
%     Nseq = Nt/Nb_perseq;
%    
%     %%
%     x_s = x;
%     x_s = reshape(x, [Ni*Nb_perseq Nseq]);
%     o_l = randperm(Nseq);
%     x_s = reshape( x_s(:, o_l), [Ni Nb_perseq*Nseq]);
%     C1 = corr_wdelay_l(x_s, Tcorr)
%     %%
%     C2 = TRBM_corr_wdelay_l(x, Nb_perseq, Tcorr)
%     
%     %%
%     figure; hold on
%     plot(C0, C1, 'k.');
%     plot(C0, C2, 'r.');
%     plot_identity_line;
%     
end