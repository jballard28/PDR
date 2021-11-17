function     [ alphahat,beta,alpha ] = simulate_sumstats_MVNmix( RR,nn,Sigma,p_causalsnp,h2,Rsqrt)
%SIMULATE_CORRELATED_PHENOTYPES_MULTI generates sumstats for three
%phenotypes with multiple shared genetic components
%
%   Input arguments: RR: LD matrix (use speye(m) for sims with no LD).
%   nn: sample size for each trait. Sigma: cell array of covariance
%   matrices, length no_cpts (each no_traits x no_traits)
%   h2g: heritability of each trait, length no_traits.
%   p_causalsnp: 1xno_cpts vector, entries sum to <=1
%   Rsqrt: optionally, provide a matrix square root of R to
%   avoid recomputing it (use speye(m) for sims with no LD).
%
%   Output arguments:
% alphahat: estimated marginal effect sizes. beta: true causal effect
% sizes. alpha: true marginal effect sizes.

no_snps=size(RR,1);
no_cpts=length(Sigma);
no_traits=length(nn);

if isscalar(p_causalsnp)
    p_causalsnp=p_causalsnp*ones(1,no_cpts);
end

if sum(p_causalsnp)>1
    error('causal proportions should sum to at most 1')
end

% For each SNP, pick a random integer corresponding to the cpt it comes
% from; no_cpts+1 -> no causal effect
choose_cpt=randsample(no_cpts,no_snps,true,p_causalsnp);

beta=zeros(no_snps,no_traits);
for ii=1:no_cpts
    beta(choose_cpt==ii,:)=mvnrnd(zeros(1,no_traits),Sigma{ii},sum(choose_cpt==ii));
end

if exist('h2')
    beta=beta*diag(sqrt(h2./sum(beta.*(RR*beta))));
end

alpha=RR*beta;

if ~exist('Rsqrt')
    Rsqrt=sqrtm(RR);
end

sumstat_noise=Rsqrt*randn(no_snps,no_traits)*diag(1./sqrt(nn));

alphahat=alpha+sumstat_noise;

end

