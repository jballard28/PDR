function output = expected_observed_r2(est,data,traitnames,alpha,replication_file,traitidx)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   replication_files: cellarray of paths to replication data file
%   traitidx: default 1

est.traits=traitnames;
alpha = alpha ./ sqrt(diag(est.cov))';

if ~isempty(replication_file)
    [rsids,z] = import_sumstat_files(replication_file);
    [~,i,j]=intersect(rsids,data.snps);
    zscores = data.z;
    ntraits=size(zscores,2);
    normalizer=mean(z(i,1).*data.z(j,traitidx));
    z = z / normalizer;
else
    zscores = data.z;
    ntraits=size(zscores,2);
    j = (1:data.noSNPs)';
end

clear *slope*
nblocks = 100;
blocksize = ceil(size(zscores,1)/nblocks);
whichBlock = ceil((1:size(zscores,1))/blocksize);

[beta_MTAG_pm,~] = MTAG_pm(zscores,data.cov,data.sigmaEps);


edges = [0, .8, .95, .99, .998, 1];

clear expected* observed*

observed_pdr_varexp = [];
observed_pdr_varexp_se = [];
observed_mtag_varexp = [];
observed_mtag_varexp_se = [];

for k1=1:ntraits
    bins = discretize(zscores(j,k1).^2,quantile(zscores(j,k1).^2,edges));
    for k2=1:length(edges)-1
        % which SNPs
        incl = bins == k2;
        incl_j = j(incl);
        
        % PDR expected
        [~,~, expected_jk] = jackknife_regression(ones(size(incl_j)), alpha(incl_j,1).^2,whichBlock(incl_j));
        expected(k1,k2) = mean(log10(expected_jk));
        expected_se(k1,k2) = std(log10(expected_jk))* sqrt(nblocks);
        
        % MTAG expected
        [~,~, mtag_expected_jk] = jackknife_regression(ones(size(incl_j)), beta_MTAG_pm(incl_j,1).^2,whichBlock(incl_j));
        mtag_expected(k1,k2) = mean(log10(mtag_expected_jk));
        mtag_expected_se(k1,k2) = std(log10(mtag_expected_jk))* sqrt(nblocks);
        
        
        if ~isempty(replication_file)
            incl_i = i(incl);
            
            % PDR observed
            [~,~, observed_jk] = jackknife_regression(ones(size(incl_j)), (alpha(incl_j,1) .* z(incl_i,1)),whichBlock(incl_j));
            observed_jk = observed_jk.^2 ./ expected_jk;
            observed(k1,k2) = mean(log10(observed_jk));
            observed_se(k1,k2) = std(log10(observed_jk)) * sqrt(nblocks);
            
            % MTAG observed
            [~,~, observed_jk] = jackknife_regression(ones(size(incl_j)), (beta_MTAG_pm(incl_j,1) .* z(incl_i,1)),whichBlock(incl_j));
            observed_jk = observed_jk.^2 ./ mtag_expected_jk;
            mtag_observed(k1,k2) = mean(log10(observed_jk));
            mtag_observed_se(k1,k2) = std(log10(observed_jk)) * sqrt(nblocks);
            
            observed_pdr_varexp = observed;
            observed_pdr_varexp_se = observed_se;
            observed_mtag_varexp = mtag_observed;
            observed_mtag_varexp_se = mtag_observed_se;
        end
        
        expected_pdr_varexp = expected;
        expected_pdr_varexp_se = expected_se;
        expected_mtag_varexp = mtag_expected;
        expected_mtag_varexp_se = mtag_expected_se;
        
    end
end

% Expected vs observed for set of all SNPs
incl_j = j;
[~,~, expected_jk] = jackknife_regression(ones(size(incl_j)), alpha(incl_j,1).^2,whichBlock(incl_j));

% MTAG expected
[~,~, mtag_expected_jk] = jackknife_regression(ones(size(incl_j)), beta_MTAG_pm(incl_j,1).^2,whichBlock(incl_j));

observed_jk = [];
mtag_observed_jk = [];
if ~isempty(replication_file)
    incl_i = i;
    
    % PDR observed
    [~,~, observed_jk] = jackknife_regression(ones(size(incl_j)), (alpha(incl_j,1) .* z(incl_i,1)),whichBlock(incl_j));
    observed_jk = observed_jk.^2 ./ expected_jk;
    
    % MTAG observed
    [~,~, mtag_observed_jk] = jackknife_regression(ones(size(incl_j)), (beta_MTAG_pm(incl_j,1) .* z(incl_i,1)),whichBlock(incl_j));
    mtag_observed_jk = mtag_observed_jk.^2 ./ mtag_expected_jk;
    
end


data.runLDscore;
noBlocks = 100;
intercept_jk = reshape(data.ldsc.intercept_jk(1,1,:),noBlocks,1);
sumstats_expected_jk_allSNPs = 1./(1 + intercept_jk);
sumstats_expected = mean(sumstats_expected_jk_allSNPs);
sumstats_expected_se = std(sumstats_expected_jk_allSNPs) * sqrt(noBlocks);

mtag_expected=mean(mtag_expected_jk);
mtag_expected_se=std(mtag_expected_jk)*sqrt(noBlocks);
expected=mean(expected_jk);
expected_se=std(expected_jk)*sqrt(noBlocks);

mtag_observed = [];
mtag_observed_se = [];
observed = [];
observed_se = [];
if ~isempty(replication_file)
    mtag_observed=mean(mtag_observed_jk);
    mtag_observed_se=std(mtag_observed_jk)*sqrt(noBlocks);
    observed=mean(observed_jk);
    observed_se=std(observed_jk)*sqrt(noBlocks);
end

output = struct('pdr_observed',observed,'pdr_observed_se',observed_se,...
    'pdr_expected',expected,'pdr_expected_se',expected_se,...
    'mtag_expected',mtag_expected,'mtag_expected_se',mtag_expected_se,...
    'mtag_observed',mtag_observed,'mtag_observed_se',mtag_observed_se,...
    'sumstats_expected',sumstats_expected,'sumstats_expected_se',sumstats_expected_se,...
    'expected_pdr_varexp',expected_pdr_varexp,'expected_pdr_varexp_se',expected_pdr_varexp_se,...
    'expected_mtag_varexp',expected_mtag_varexp,'expected_mtag_varexp_se',expected_mtag_varexp_se,...
    'observed_pdr_varexp',observed_pdr_varexp,'observed_pdr_varexp_se',observed_pdr_varexp_se,...
    'observed_mtag_varexp',observed_mtag_varexp,'observed_mtag_varexp_se',observed_mtag_varexp_se,...
    'pdr_expected_jk',expected_jk,'pdr_observed_jk',observed_jk,...
    'mtag_expected_jk',mtag_expected_jk,'mtag_observed_jk',mtag_observed_jk);

end

