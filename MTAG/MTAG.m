function [beta_MTAG,MTAG_weights] = MTAG(Z,Rg,sigmaEps)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[n,k]=size(Z);

for j=1:k
    rj=Rg(:,j);
    MTAG_weights(j,:) = (rj' / (Rg-rj*rj'+sigmaEps)) / ...
        ((rj' / (Rg-rj*rj'+sigmaEps)) * rj);
end

beta_MTAG = (MTAG_weights*Z')';

end

