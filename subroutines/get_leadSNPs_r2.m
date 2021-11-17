function [leadSNPs,h2GWAS] = get_leadSNPs_r2(chisq,RRb,sig_thresh)
%get_leadSNPs_r2 performs greedy LD pruning, iteratively discarding all SNPs in LD with lead SNPs. 
%   chisq: pruning proceeds from largest to smallest value on this list.
%   RRb: sparse matrix of which SNPs are in LD with each other (e.g.
%   r^2>0.01). sig_thresh: chisq threshold beyond which to stop reporting
%   results.

incl=chisq>sig_thresh;
RRb=RRb(incl,incl);
chisq=chisq(incl);

[chisq,ix]=sort(chisq,'descend');
RRb=RRb(ix,ix);

ii=1;
incl=true(size(chisq));
leadSNPs=false(size(chisq));
while any(incl)
    ii=ii+1;
    imax=find(incl,1,'first');
    incl(RRb(:,imax))=false;
    leadSNPs(imax)=true;
end

h2GWAS=sum(chisq(leadSNPs));
end

