function [c] = sparsecell(ii,jj,vv,n,m,default)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
c=repmat({default},n,m);
c(ii+(jj-1)*n)=vv;
end

