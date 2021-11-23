function make_scatterplot(est,data,alpha,posterior_scalars,cptnames,traits,subplot_traits,included_cpts,LDSNPs,RRb)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% LD pruning
[~,i,j]=intersect(data.snps,LDSNPs);
SNPorder = zeros(length(i),est.noTraits);
for tt=1:est.noTraits
    [~, SNPorder(:,tt)] = sort(alpha(i,tt).^2);
end
SNPorder = max(SNPorder')';
leadSNPs = get_leadSNPs_r2(SNPorder,RRb(j,j),0.99 * length(SNPorder));
incl = i(leadSNPs);


% Rescaling factor to help assign SNPs to components
clear rescale h2Cpt
for i=1:est.noCpts
    rescale(i)=trace(est.cpts((i)).cov)/sum(posterior_scalars(:,(i)));
end
scalars = posterior_scalars(incl, :) .* rescale;
pm=alpha(incl,:);

% Assign each SNP to its component with the largest normalized scalar; any
% components not included in included_cpts are merged into 'other'
[~,cptAssignment] = max(scalars');
reassignment = ones(1,est.noCpts) * (length(included_cpts) + 1);
reassignment(included_cpts) = 1:length(included_cpts);
assignment = reassignment(cptAssignment);

% Colors to plot
colors=colororder;
colors = colors(1:length(included_cpts),:);
colors(end+1,:) = [.5 .5 .5];


%% scatterplots
axlim=quantile(abs(pm(:)),.999);

for k=1:size(subplot_traits,1)
    figure;
    hold on
    scatter(pm(:,subplot_traits(k,1)),...
        pm(:,subplot_traits(k,2)),3,'filled','CData',colors(assignment,:),'Marker','o');
    for ii = 1:length(included_cpts)
        h(ii)=scatter(nan,nan,1,colors(ii,:),'filled');
    end
    legend(h, cptnames)
    xlim(axlim*[-1 1]);
    ylim(axlim*[-1 1]);
    set(gca,'XTick',0,'XTickLabel',[])
    set(gca,'YTick',0,'YTickLabel',[])
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');
    box off
    title(sprintf('%s vs %s\n(r_g=%.2f)',traits{subplot_traits(k,2)},...
        traits{subplot_traits(k,1)},data.cov(subplot_traits(k,1),subplot_traits(k,2))));
end



end

