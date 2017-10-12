function [M, dim_l] = shiftdim_nosqueeze(M, n, dim_l )

if nargin < 3
    dim_l = size(M);
end
%         dim_l = size(M);
D = length(dim_l);
if n<1
    n = mod(n-1,D)+1;
end

if n ~= D
    o_l = mod((0:(D-1)) + n,D) +1; % equivalent but faster than [(n+1:D), 1:n];
    M = permute(M, o_l);
    dim_l = dim_l(o_l);
    
end

% if 0
%     %% test
%     D = 3;
%     n = 2;
%     tic
%     for k = 1:1000
%         o_l = [(n+1:D), 1:n];
%     end
%     toc
%     
%     tic
%     for k = 1:1000
%         o_l = mod((0:(D-1)) + n,D) +1;
%     end
%     toc
% end
end

