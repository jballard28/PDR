function [beta_MTAG,MTAG_weights] = MTAG_pm(Z,Rg,sigmaEps)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[n,k]=size(Z);

MTAG_weights = ((inv(Rg)+inv(sigmaEps))\inv(sigmaEps));

beta_MTAG = (MTAG_weights*Z')';

end

