function plot_varexp(r2_output,traitnames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ntraits = length(traitnames);

pdr_expected_jk = r2_output.pdr_expected_jk;
pdr_observed_jk = r2_output.pdr_observed_jk;
mtag_observed_jk = r2_output.mtag_observed_jk;
mtag_expected_jk = r2_output.mtag_expected_jk;
expected_pdr_varexp = r2_output.expected_pdr_varexp;
expected_pdr_varexp_se = r2_output.expected_pdr_varexp_se;
expected_mtag_varexp = r2_output.expected_mtag_varexp;
expected_mtag_varexp_se = r2_output.expected_mtag_varexp_se;
observed_pdr_varexp = r2_output.observed_pdr_varexp;
observed_pdr_varexp_se = r2_output.observed_pdr_varexp_se;
observed_mtag_varexp = r2_output.observed_mtag_varexp;
observed_mtag_varexp_se = r2_output.observed_mtag_varexp_se;

nblocks=100;

% PDR expected_pdr_varexp
expected_allSNPs = mean(log10(pdr_expected_jk));
expected_se_allSNPs = std(log10(pdr_expected_jk))* sqrt(nblocks);

% PDR observed_pdr_varexp
observed_allSNPs = mean(log10(pdr_observed_jk));
observed_se_allSNPs = std(log10(pdr_observed_jk)) * sqrt(nblocks);

% MTAG expected_pdr_varexp
mtag_expected_allSNPs = mean(log10(mtag_expected_jk));
mtag_expected_se_allSNPs = std(log10(mtag_expected_jk))* sqrt(nblocks);

% MTAG observed_pdr_varexp
mtag_observed_allSNPs = mean(log10(mtag_observed_jk));
mtag_observed_se_allSNPs = std(log10(mtag_observed_jk)) * sqrt(nblocks);


figure;subplot(1,2,1);hold on
plot([-1,2], [-1,2], 'black');
for k1 = 1:ntraits
    h(k1) = errorbar(expected_pdr_varexp(k1,:),observed_pdr_varexp(k1,:),...
        observed_pdr_varexp_se(k1,:),observed_pdr_varexp_se(k1,:),...
        expected_pdr_varexp_se(k1,:),expected_pdr_varexp_se(k1,:),...
        '.','MarkerSize',12,'CapSize',0);
end
h(k1+1) = errorbar(expected_allSNPs,observed_allSNPs,...
    observed_se_allSNPs,observed_se_allSNPs,...
    expected_se_allSNPs,expected_se_allSNPs,...
    '.','MarkerSize',12,'CapSize',0,'Color','black');
legend_names = cellfun(@(s){[s,' bins']},traitnames);
legend_names(end+1)={'All SNPs'};
legend(h([end,1:end-1]), legend_names([end,1:end-1]));legend boxoff
xlabel('Expected variance'),ylabel('Observed variance')
tickValues = log10(10.^(-1:2).*[1 2 5]');
set(gca,'XTick',tickValues(:),'XTickLabel',{'0.1','','','1','','','10','','','100','',''},...
    'XLim',[-1-log10(2),2+log10(2)],'YLim',[-1-log10(2),2+log10(2)]);
set(gca,'YTick',tickValues(:),'YTickLabel',{'0.1','','','1','','','10','','','100','',''});
set(gca,'TickLength',[.02 .02])
ax = gca;
title('PDR per-SNP variance explained')

subplot(1,2,2);hold on
plot([min(expected_pdr_varexp(:)),max(expected_pdr_varexp(:))], [min(expected_pdr_varexp(:)), max(expected_pdr_varexp(:))], 'black');
for k1 = 1:ntraits
    h(k1) = errorbar(expected_mtag_varexp(k1,:),observed_mtag_varexp(k1,:),...
        observed_mtag_varexp_se(k1,:),observed_mtag_varexp_se(k1,:),...
        expected_mtag_varexp_se(k1,:),expected_mtag_varexp_se(k1,:),...
        '.','MarkerSize',12,'CapSize',0);
end
h(k1+1) = errorbar(mtag_expected_allSNPs,mtag_observed_allSNPs,...
    mtag_observed_se_allSNPs,mtag_observed_se_allSNPs,...
    mtag_expected_se_allSNPs,mtag_expected_se_allSNPs,...
    '.','MarkerSize',12,'CapSize',0,'Color','black');

xlabel('Expected variance')%ylabel('Observed variance')
set(gca,'XLim',ax.XLim,'YLim',ax.YLim,'XTick',ax.XTick,'YTick',ax.YTick,...
    'XTickLabel',ax.XTickLabel,'YTickLabel',ax.YTickLabel,'TickLength',ax.TickLength)
title('MTAG per-SNP variance explained')

end

