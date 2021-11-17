%% get prediction with posterior variance (posterior_scalars)

% (1) load the model
% load(modeldir)

% (2) call predict (optional)
[alpha, xpm, likelihood, posterior_scalars] = est.predict(data145,'minWeight',1/(1e6*data145.noSNPs));

% (3) save output of predict (optional)
save('./gene_enrichment/data/realtraits_060821_T2D-BMI-FG-A1C_2fr_predresults','alpha', 'xpm', 'likelihood', 'posterior_scalars')

% (4) load prediction results
load('./gene_enrichment/data/realtraits_060821_T2D-BMI-FG-A1C_2fr_predresults')

trait_idx = 1;

%% get SNPs that surpass a given threshold for all components

percentile=95;

% Merging all posterior scalars across components to find universal threshold
for cc = 1:est.noCpts
    mult_element = est.cpts(cc).S(trait_idx,trait_idx);
    scaled_ps(:,cc) = posterior_scalars(:,cc)*mult_element;
end

all_ps = scaled_ps(:);

ps_thresh = prctile(all_ps,percentile);

% Iterating over all cpts and taking SNPs that surpass thresh
for cc = 1:est.noCpts
    
    ps = scaled_ps(:,cc);
    
    cptsnps = find(ps > ps_thresh);
    cptrs = data145.snps(cptsnps);
    
    % appending rs to the beginning of each id
    cptrs = string(arrayfun(@(x)join(['rs',string(x)],''),cptrs,'UniformOutput',false));
    
    % exporting to csv
    writematrix(cptrs,join(['./gene_enrichment/data/coronary_3fr_082621_cpt',num2str(cc),'.csv']));
    
    allcptsnps{cc} = cptrs;
    
end
