function [ corr_l ] = corr_wdelay_l( x, Tcorr )
%CORR_WDELAY_L 
% cross correlation
% corr_l(:,k) is the list of cross-correlations with delay (k-1)
% can be inverted by cov_wdelay2mat

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

corr_l = [];
if ~verLessThan('matlab', '8.3')
    for t = 0:(Tcorr-1)
%                 corr_l(:,t+1) = triu_l_withdiag(corr(x', circshift(x,-t,2)'))';
        m = corr(x', circshift(x,-t,2)');
        corr_l(:,t+1) = m(:);
            step_notification(t, 1, 15);
    end
else
    for t = 0:(Tcorr-1)
%          corr_l(:,t+1) = triu_l_withdiag(corr(x', circshift(x,-t,2)'))';
        m = corr(x', circshift_asOLD(x,-t,2)');
        corr_l(:,t+1) = m(:);
%             step_notification(t, 1, 15);
    end
end
%% ensure exactly 1 for t = 0 and  i = j 
ke =  eye(size(x,1)); %triu_l_withdiag(eye(size(x,1)));
corr_l(ke>0,1)=1;



% if 0
%     %% test
%     x = rand(2,10);
%     x_u = [];
%     Tcorr = 3
%     for t = 0:(Tcorr-1)
%         x_u = [x_u; circshift(x,-t,2)]; 
%     end
%     x_u;
%     c1 = corr(x_u')
%     %%
%     c2 = corr_wdelay_l( x, Tcorr)
%     %%
%     c2 == 1
%     
%   
% end
end
