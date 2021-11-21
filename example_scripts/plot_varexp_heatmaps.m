%% get prediction with posterior variance (posterior_scalars)

load('fitmodel.mat')

% Calling predict function to get posterior scalars
[alpha, MAP_scalars, posterior_scalars] = est.predict(data,'minWeight',1e-6);

alpha = alpha ./ sqrt(diag(est.cov))';

traits = {'T2D','BMI','TG'};
cptnames = {'Cpt 1','Cpt 2'};
included_cpts = 4:5;

% Download this file from dropbox.com/sh/mclm1urkxs8ga80/AADDDABQYeGtyQxmom2raMkva?dl=0
% wget https://www.dropbox.com/sh/mclm1urkxs8ga80/AABLDZRREAkGj5A1D3x8z_FOa/1kg_LD.HM3.window1cm.noblocks.mat?dl=0

% NOTE: modify this path
load('/Users/jballard/Documents/data/1kg_LD.HM3.window1cm.noblocks.mat','LDSNPs','RRb');

%% Which SNPs to include

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
snps=data.snps(incl);

% Assign each SNP to its component with the largest normalized scalar; any
% components not included in included_cpts are merged into 'other'
[maxsc,cptAssignment] = max(scalars');
reassignment = ones(1,est.noCpts) * (length(included_cpts) + 1);
reassignment(included_cpts) = 1:length(included_cpts);
assignment = reassignment(cptAssignment);


%% h2 plot
% Colors to plot
colors=colororder;
colors = colors(1:length(included_cpts),:);
colors(end+1,:) = [.5 .5 .5];

figure;
est.h2plot('traitNames',traits,'cptNames',cptnames,'whichCpts',included_cpts);
title('Variance explained by each component')
