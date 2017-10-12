function [ cov_mat ] = cov_wdelay2mat( cov_l, T)
% transforms cov_l into a matrix of size
% (N, T, N, T) of correlations between (x_(i,tau+t), x_(j, tau + t'))

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

[N,B] = size(cov_l);
N = sqrt(N);

if ~exist('T','var')
    T=B;
end

cov_mat = zeros(N,N,B,B);
for ti = 1:T
    for tj = ti:T
        if (tj-ti+1)<=B
            cov_mat(:,:,tj,ti) = reshape(cov_l(:,tj-ti+1),[N N]);
            cov_mat(:,:,ti,tj) =  cov_mat(:,:,tj,ti)';
        end
    end
end
cov_mat = permute(cov_mat,[1 3 2 4]);
cov_mat = reshape(cov_mat,[N*T, N*T]);




if 0
    %%
    N = 4;
    T = 200;
    %     x = rand(N,T);
    x = rand(1,T);
    x(2,2:end) = 2*x(1:(end-1));
    
    x = [x(:,2:end);  -1.5*x(:,1:(end-1))];
    
    D = 3
    %%
    b = TRBM_unfold_time(x,Inf,D);
    cov_mat = cov( b');
    
    cov_l = cov_wdelay_l(x, D);
    cov_mat = cov_wdelay2mat( cov_l)
    
    isequal(Cc, cov_mat)
    %%
    T = 4
    cov_mat = cov_wdelay2mat( cov_l, T)

    %%
    figure; imagesc(cov_mat == 0)
    %%
    figure; imagesc(Cc)
    
    %%
    figure; hold on
    plot(unique(cov_mat(:)));
    
    plot(unique(Cc(:)));
    %%
    figure; plot(Cc(:), cov_mat(:),'.');
    plot_identity_line;
end

