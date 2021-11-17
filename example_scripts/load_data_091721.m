%% Loading data

% clustertraits: cellarray with names of sumstats files
% clustername: name of trait cluster
% tnames: cellarray with trait names in same order as clustertraits

% Load data
l2fp='/Users/jballard/Documents/LDSC/tutorial/baseline_unzipped/baselineLD.';
weightsfp='/Users/jballard/Documents/LDSC/tutorial/1000G_Phase3_weights_hm3_no_MHC_unzipped/weights.hm3_noMHC.';
fp='/Users/jballard/Documents/data/';

data=DATA(cellfun(@(s){[fp,s,'.sumstats']},clustertraits),'LDscoreDir',l2fp,'WeightsDir',weightsfp);

est_data = copy(data);