function  LDSCout = CLDSC( data,varargin)
%CSLD4M (cross-trait stratified LD 4th moments regression).
%
%   REQUIRED INPUT ARGUMENTS:chisq1,chisq2:  Mx1 vectors of chi^2 statistics;
%   zprod: Mx1 vector which is product of Z scores;
%   data.l2: MtotxP matrix of LD scores (LD second moments), where P is #data.annotations
%   and Mtot is the number of reference SNPs;
%   l4: MtotxP matrix of LD 4th moments;
%   Optionally, prepend an all-ones vector in first column of data.l2 and/or l4,
%   which will cause an intercept to be used in the regression.
%   data.annot: MtotxP matrix of data.annotation values;

%   OPTIONAL INPUT ARGUMENTS as name-value pairs:
%   RegressionIndices: Mx1 vector of indices corresponding to the reference SNPs that
%   will be used in the regression (regression SNPs must be a subset of
%   reference SNPs);
%   data.l2Weights: Mx1 vector of weights for LDSC;
%   NoJackknifeBlocks: number of jackknife blocks (default 100);
%   OutputSubspace: which data.annotations to report estimates for. This can be
%   either a list of P' indices <=P, or a PxP' projection matrix whose columns
%   correspond to linear combinations of input data.annotations.
%
%
%   OUTPUT ARGUMENTS: all outputs are matrices of size NoJacknifeBlocks x
%   P', ie each row = a leave-one-block-out estimate and each column = an
%   data.annotation.
%   rg: genetic correlation
%   h21: "heritability" for trait 1, if input summary stats are in
%   standardized units
%   h22: "heritability" for trait 2, if input summary stats are in
%   standardized units

mm_regression=length(data.z);

p=inputParser;

addRequired(p,'data',@(x)isa(x,'DATA') && isscalar(x))
addOptional(p,'NoJackknifeBlocks',100,@(x)isscalar(x) && all(mod(x,1)==0) && all(x<=mm_regression) && all(x>1))
addOptional(p,'OutputSubspace',[],@(x)checkProjectionmatrix(x) || checkIndexlist(x))
addOptional(p,'BaseCol',1,@isscalar)
addOptional(p,'noSNPsAnnot',[],@true)

parse(p,data,varargin{:});

no_blocks=p.Results.NoJackknifeBlocks;
noSNPsAnnot = p.Results.noSNPsAnnot;

if isempty(data.weights)
    l2weights=1;
else
    l2weights=1./data.weights;% 2nd-moment weights matrix
end
WW2_mat=diag(sparse(l2weights));

clear varargin

[mm,pp]=size(data.l2);

% add intercept column to l2
l2 = [ones(mm,1), data.l2];
pp=pp+1;

if mm~=mm_regression
    error('Wrong number of snps')
end

blocksize=floor(mm/no_blocks);

ell_ell=zeros(pp,pp,no_blocks);

tp = data.noTraits * (data.noTraits+1) / 2;
ell_zz=zeros(pp,tp,no_blocks);

zprod = zeros(mm,tp);
counter=1;
for i=1:data.noTraits
    for j=1:i
        zprod(:,counter) = data.z(:,i).*data.z(:,j);
        counter = counter+1;
    end
end

% Compute various moments for each jackknife block
for jk=1:no_blocks
    ind=(jk-1)*blocksize+1:jk*blocksize;
    ell_ell(:,:,jk)=l2(ind,:)'*WW2_mat(ind,ind)*l2(ind,:)/sum(data.weights(ind));
    ell_zz(:,:,jk)=l2(ind,:)'*WW2_mat(ind,ind)*zprod(ind,:)/sum(data.weights(ind));%
end

% Combine jackknife estimates to obtain leave-one-out regression
% coefficients
tau=zeros(pp,tp,no_blocks);

for jk=1:no_blocks
    ind=[1:jk-1,jk+1:no_blocks];
    tau(:,:,jk) = mean(ell_ell(:,:,ind),3)\mean(ell_zz(:,:,ind),3);% data.l2var^-1 * data.l2a2
    if ~isempty(data.noSNPsAnnot)
        h2total(:,jk) = data.noSNPsAnnot * tau(2:end,:,jk) ;
    end
end

intercept = tau(1,:,:);
LDSCout.intercept = triuind(mean(intercept,3),1);
LDSCout.intercept_se = triuind(std(intercept,[],3)*sqrt(no_blocks),1);
LDSCout.intercept_jk = intercept;

tau = tau(2:end,:,:);
LDSCout.tau=mean(tau,3);
LDSCout.tau_se = std(tau,[],3)*sqrt(no_blocks);
LDSCout.tau_jk = tau;

perSNPh2 = tau(p.Results.BaseCol,:,:);
LDSCout.perSNPh2=triuind(mean(perSNPh2,3),1);
LDSCout.perSNPh2_se=triuind(std(perSNPh2,[],3)*sqrt(no_blocks),1);

if ~isempty(data.noSNPsAnnot)
    Rg = triuind(mean(h2total,2));
    D = diag(Rg);
    Rg = Rg + Rg' - diag(D);    
    LDSCout.rg = diag(sqrt(1./D)) * Rg * diag(sqrt(1./D));
    
end
%LDSCout.rg=corrcov(LDSCout.perSNPh2);


end




