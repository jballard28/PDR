function pval_qq_plot(null_pvals)

% QQ plot for pvals
fontsize=16;
markersize=10;

unif_quant=-log10(rand(size(null_pvals,2),1));
qq1=qqplot(unif_quant,-log10(null_pvals(1,:)));
qq1(1).MarkerEdgeColor='r';
hold on;
qq2=qqplot(unif_quant,-log10(null_pvals(2,:)));
qq2(1).MarkerEdgeColor='b';
delete(qq1(2))
delete(qq1(3))
delete(qq2(3))
delete(qq2(2))

qq1(1).Marker = '.';
qq2(1).Marker = '.';

qq1(1).MarkerSize=markersize;
qq2(1).MarkerSize=markersize;

xx=[0 max(unif_quant)];
hold on;
p1=plot(xx,xx,'k--');
legend([qq1(1) qq2(1) p1],{'null','non-null','y=x'},'Location','northwest','FontSize',12);

xlabel('-log10(Uniform Quantiles)','FontSize',fontsize);
ylabel('-log10(1 pleio cpt vs. 2 pleio cpts pval)','FontSize',fontsize);

end

