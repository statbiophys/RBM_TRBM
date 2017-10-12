function L = triu_l( m )
%TRIU_L list of the elements above the diagonal of m

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

if length(size(m))==2
    if diff(size(m))
        error('m must be squarre');
    end
    
    L = m(triu(ones(size(m))>0,1));
else
    s_l = size(m);
    if s_l(1) ~= s_l(2)
        error('m must be squarre in first 2 dimensions');
    end
    ke = triu(ones(s_l(1))>0,1);
    m = reshape(m,[s_l(1)*s_l(2), s_l(3:end)]);
    L = m(ke(:),:);
    L = L(:);
end

% if 0
%     %%
%     m = rand(3,3,2)%>0.5
%     triu_l(m)
% end
end

