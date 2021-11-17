function [ beta,exitflag ] = nonnegative_regression( X,y,W,init,ACon,VCon )
%Nonnegative linear regression of y on X with weights W
%   init: initial value for optimization (should satisfy constraints).
%   ACon, VCon: additional equality constraints that should be satisfied,
%   ie ACon*beta==VCon

[n,p]=size(X);
if ~exist('W')
    W=ones(n,1);
end
if isvector(W)
   W=diag(sparse(W));
end
if nargin<4
    init=zeros(size(X,2),1);
elseif isempty(init)
        init=zeros(size(X,2),1);
end
if nargin<6
    ACon=[];
    VCon=[];
end
if p>n
    warning('This method of solving nonnegative least squares is inefficient for underdetermined systems')
end

% alpha=X'*W*y/n;
% Sigma=X'*W*X/n;% pxp

alpha=y;
Sigma=X;% pxp

% Sigma(Sigma == Inf | isnan(Sigma)) = 1e10;
% Sigma(Sigma == -Inf) = -1e10;
% Sigma
%options = optimset('MaxIter',1000000);
%[beta,~,~,exitflag]=lsqnonneg(Sigma,alpha,options);

%options = optimoptions('lsqlin','Display','off');
[beta,~,~,exitflag]=lsqlin(Sigma,alpha,...
    -eye(p),zeros(p,1),...
    ACon,VCon,...
    [],[],init,optimoptions('lsqlin','Display','off'));% Regress alpha on Sigma with constraint -eye(p)*beta<=zeros(p,1)
% sum((Sigma*beta-alpha).^2)
%beta=(X'*diag(sparse(W))*X)\(X'*diag(sparse(W))*y);

% Setting betas to zero that when zero decrease the objective function
% sum of squared residuals
cur_obfun_val=dot(alpha-Sigma*beta,alpha-Sigma*beta);
set0_idcs=[];
for i=1:length(beta)
   temp_beta=beta;
   temp_beta(i)=0;
   temp_beta=temp_beta*sum(beta)/sum(temp_beta);
   temp_obfun_val=dot(alpha-Sigma*temp_beta, alpha-Sigma*temp_beta);
   if temp_obfun_val <= cur_obfun_val
       set0_idcs = [set0_idcs i];
   end
end
% set0_idcs
% beta(set0_idcs)
beta(set0_idcs) = 0;

end

