function [ mEv_l] = TRBM_cyclic_h2mEv( M, h )
% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

Tmem = size(M.w,3);

tm = nmodeproduct(M.w, h', 1);
t = tm(:,:,1);
if ~verLessThan('matlab', '8.3')
    for d = 1:(Tmem-1)
        t = t + circshift(tm(:, :, d+1),-d, 1);
    end
else
    for d = 1:(Tmem-1)
        t = t + circshift_asOLD(tm(:, :, d+1),-d, 1);
    end
end

mEv_l =  bsxfun(@plus, squeeze(t)', M.a' );

% if 0
%     %% test
%     Nt = 10;
%     Nj = 3;
%     Ni = 4;
%     Tmem = 8;
%     h = rand(Nj, Nt);
%     M.w = rand(Nj, Ni, Tmem);
%     tm =  nmodeproduct(M.w, h', 1);
%     
%     
%     %% much slower
%     z = permute(tm,[2 3 1]);
%     d = size(z)
%     z = reshape(z, [d(1)*d(2), d(3)])
%     [a,~,c] = find(z)
%     o = repmat(0:(d(2)-1),[d(1), 1]);
%     t = bsxfun(@plus, o(:)  ,1:d(3))
%     
%     t2 = full(sparse(a, t(:), c));
%     t2 = reshape(t2, [d(1), d(2), d(3)+Tmem-1])
%     
%     t2 = full(t2)
%     t2 = find(tm);
%     
%     %%
%     isequal(t,t2)
%     figure; hist(t(:) - t2(:));
%     figure; plot(t(:), t2(:),'.');
end