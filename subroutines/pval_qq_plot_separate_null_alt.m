function pval_qq_plot_separate_null_alt(null_pvals)

% QQ plot for pvals
fontsize=16;
markersize=10;

unif_quant=-log10(rand(size(null_pvals,2),1));

figure;
set(gcf,'Position',[458 492 1080 454]);
subplot(1,2,1);
qq1=qqplot(unif_quant,-log10(null_pvals(1,:)));
qq1(1).MarkerEdgeColor='r';
delete(qq1(2))
delete(qq1(3))
qq1(1).Marker = '.';
qq1(1).MarkerSize=markersize;

xx=[0 max(unif_quant)];
hold on;
p1=plot(xx,xx,'k--');
legend([qq1(1) p1],{'null','y=x'},'Location','northwest','FontSize',12);

xlabel('-log10(Uniform Quantiles)','FontSize',fontsize);
ylabel('-log10(1 pleio cpt vs. 2 pleio cpts pval)','FontSize',fontsize);

subplot(1,2,2);
qq2=qqplot(unif_quant,-log10(null_pvals(2,:)));
qq2(1).MarkerEdgeColor='b';
delete(qq2(3))
delete(qq2(2))
qq2(1).Marker = '.';
qq2(1).MarkerSize=markersize;

xx=[0 max(unif_quant)];
hold on;
p1=plot(xx,xx,'k--');
legend([qq2(1) p1],{'non-null','y=x'},'Location','northwest','FontSize',12);

xlabel('-log10(Uniform Quantiles)','FontSize',fontsize);
ylabel('-log10(1 pleio cpt vs. 2 pleio cpts pval)','FontSize',fontsize);

end

