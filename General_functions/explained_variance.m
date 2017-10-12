function [ ev, Numerator, Denominator, Denominator_ste  ] = explained_variance( data_l, model_l )
%EXPLAINED_VARIANCE
% explained variance of data explained by model

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

Denominator = var(data_l);
Numerator = Denominator -  mean( (data_l - model_l).^2 );

ev = Numerator/Denominator;

if nargout>3
    l = (data_l - mean(data_l)).^2;
    Denominator_ste = std(l)/sqrt(length(l));
end
end

