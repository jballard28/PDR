function covdiff = covdiff(cov1,cov2)

covdiff = abs(cov1-cov2)./abs(cov2);
covdiff = mean(covdiff,'all');

end