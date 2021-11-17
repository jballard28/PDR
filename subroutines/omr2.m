function omr2 = omr2(z1,z2)
% Calculate 1-r^2 between two input effect size vectors

    omr2 = diag(1-corr(z1,z2));

end