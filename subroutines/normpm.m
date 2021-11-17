function [xpm, likelihood] = normpm(x,sigma,noisevar)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


likelihood = mvnpdf(x,zeros(size(x)),sigma + noisevar);

if det(sigma) < 1e-12
    sigma = sigma + 1e-12 * eye(size(sigma));
end

xpm = (inv(noisevar) + inv(sigma)) \ (noisevar\x');

end

