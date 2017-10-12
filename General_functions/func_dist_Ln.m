function v = func_dist_Ln( arg_1, arg_2, n )
%FUNC_DIST_Ln(arg_1, arg_2, n)
% L_n norm of the difference between arg_1 and arg_2

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

% if  n == 2 
%     v = sqrt(sum((arg_1(:) - arg_2(:)).^2));
% else % faster
    v = norm(arg_1(:)-arg_2(:),n);
% end
end


