function [slope_est,slope_se,slope] = jackknife_regression(X,y,whichBlock)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[nrow,ncol]=size(X);

if nargin<3
    nblocks = 100;
    blocksize = ceil(nrow/nblocks);
    whichBlock = ceil((1:nrow)/blocksize);
else
    nblocks=max(whichBlock);
end

XtX=zeros(ncol,ncol,nblocks);
XtY=zeros(ncol,nblocks);
for jk = 1:nblocks
    XtX(:,:,jk) = X(whichBlock == jk,:)'*X(whichBlock == jk,:);
    XtY(:,jk) =  X(whichBlock == jk,:)'*y(whichBlock == jk);
end
XtXsum=sum(XtX,3);
XtYsum=sum(XtY,2);

slope = zeros(ncol,nblocks);
for jk=1:nblocks
    slope(:,jk)=(XtXsum-XtX(:,:,jk)) \ (XtYsum - XtY(:,jk));
end

slope_est=mean(slope,2);
slope_se = std(slope')' * sqrt(nblocks+1);
