function y = sigmoid_approx(x)
% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

if numel(x)<2048 
    y = mex_sigmoid_fast(x); % faster for short arrays
else
    y = sigmoid(x); % faster for long arrays
end

% if 0
%     %% test
%     t1 = []; t2 = [];
%     L = ceil(2042:1:2052);
%     for l = L
%         x =  normrnd(1.5,2,1,l);
%         tic
%         for  k = 1:10000
%             y3 = mex_sigmoid_fast(x);
%         end
%         t2(end+1) = toc;
%         tic
%         for  k = 1:10000
%             y3 = sigmoid(x);
%         end
%         t1(end+1) = toc;
%         l
%     end
%     toc
%     %%
%     figure; hold on
%     plot(L, t1./L);
%     plot(L, t2./L);
%     ax =gca;
%     ax.XScale = 'log';
%     %%
%     figure; hold on
%     x = -5:0.01:5;
%     plot(x, sigmoid(x)*100);
%     plot(x, mex_sigmoid_fast(x)*100);
% end
end
