%% Loading data

clustertraits = {'/Users/jballard/Documents/data/PASS_Type_2_Diabetes.sumstats','/Users/jballard/Documents/data/PASS_BMI1.sumstats','/Users/jballard/Documents/data/PASS_Triglycerides.sumstats'};
l2fp = '/Users/jballard/Documents/LDSC/tutorial/baseline_unzipped/baselineLD.';
weightsfp = '/Users/jballard/Documents/LDSC/tutorial/1000G_Phase3_weights_hm3_no_MHC_unzipped/weights.hm3_noMHC.';

data = DATA(clustertraits, 'LDscoreDir', l2fp, 'WeightsDir', weightsfp);