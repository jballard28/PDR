function [ beta,exitflag ] = constrained_regression( X,y,W,nonneg_beta,init,ACon,VCon )
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
if nonneg_beta
    A=-eye(p);
    c=zeros(p,1);
else
    A=[];
    c=[];
end
[beta,~,~,exitflag]=lsqlin(Sigma,alpha,...
    A,c,...
    ACon,VCon,...
    [],[],init,optimoptions('lsqlin','Display','off'));% Regress alpha on Sigma with constraint -eye(p)*beta<=zeros(p,1)
% sum((Sigma*beta-alpha).^2)
%beta=(X'*diag(sparse(W))*X)\(X'*diag(sparse(W))*y);

end

