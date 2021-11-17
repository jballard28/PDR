function [grad] = gradient_function(tt,beta0,sigma,dsigma,noisevar)
%GRADIENT_FUNCTION computes partial derivative of X*beta wrt theta where X=
%phi_theta(tt) and phi_theta is a vector of N(0,sigma_theta) CFs
%   Input arguments: tt: K x no_traits vector of sampling times
%   beta0: P x 1 vector of effect sizes from previous iteration
%   theta0: 1 x Q vector of parameter estimates from previous iteration
%   delsigma: P x Q cell array of functions. Functions input a 1 x Q vector
%   of parameters (theta) and output a no_traits x no_traits matrix, the
%   partial derivative of the p-th covariance matrix wrt the q-th parameter
%   Output: d/dtheta X(theta)*beta, a K x 1 vector


f=@(t)beta0'*(cellfun(@(dS)-1/2*t*dS*t',dsigma)'.*...
    cellfun(@(S)exp(-1/2*t*(S+noisevar)*t'),sigma));
%g=@(t)beta0'*(cellfun(@(fn)exp(-1/2*t*fn(theta0)*t'),sigma)');

grad=arrayfun(@(row){f(tt(row,:))},1:size(tt,1));% use arrayfun to avoid computing unneccessary cross-terms
grad=vertcat(grad{:});


end

