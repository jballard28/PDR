# Pleiotropic Decomposition Regression (PDR)

![Method schematic](https://github.com/jballard28/PDR/blob/main/imgs/pdr_schematic.png)

# Demo
## Loading Data
You will need summary statistics that have the following columns in this exact order:
1. SNP: SNP rsID
2. A1: effect allele
3. A2: other allele
4. N: sample size
5. CHISQ: squared Z-score
6. Z: Z-score (note that it's sign is with respect to A1)

Download the summary statistics that will be input to PDR. For some example data, one can run the commands:
```
wget https://storage.googleapis.com/broad-alkesgroup-public/sumstats_formatted/PASS_Type_2_Diabetes.sumstats
wget https://storage.googleapis.com/broad-alkesgroup-public/sumstats_formatted/PASS_BMI1.sumstats
wget https://storage.googleapis.com/broad-alkesgroup-public/sumstats_formatted/PASS_Triglycerides.sumstats
```

Download LD scores and weights as follows:
```
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz
```

Decompress the LD scores and weights `.tgz` files, and then unzip each of the files within those directories. If the LD scores are contained within a directory called `baseline_unzipped`, and the weights are in a directory called `1000G_Phase3_weights_hm3_no_MHC_unzipped`, the data can be loaded and preprocessed using the following commands within MATLAB:
```
clustertraits = {'/path/to/PASS_Type_2_Diabetes.sumstats','/path/to/PASS_BMI1.sumstats','/path/to/PASS_Triglycerides.sumstats'};
l2fp = '/path/to/baseline_unzipped/baselineLD.';
weightsfp = '/path/to/1000G_Phase3_weights_hm3_no_MHC_unzipped/weights.hm3_noMHC.';

data = DATA(clustertraits, 'LDscoreDir', l2fp, 'WeightsDir', weightsfp);
```

## Fitting a PDR model
See `example_scripts/fit_model.m` for an example.

There are many parameters in this script that can be altered according to which type of model you would like to fit. For example, the number of rank-one ('factor-like') and full-rank ('generic pleiotropic') components can be altered. If more traits or more components are being fit, we recommend increasing the number of initializations (`ninit`) and also the maximum number of gradient descent steps (`gdsteps`). If the code is giving warnings that there is an insufficient number of sampling times compared to the number of parameters, we recommend reducing the `rotation_tolerance` parameter in the `initialize_fit` function within `MODEL`.

## Visualizing the components as heat maps
The `plot_heatmap` function (in `plotting`) visualizes heatmaps of the component covariance matrices that are normalized such that the total variance explained across all components for each trait is equal to 1. Below is an example of how to run this function using the output of `fit_model.m`:
```
tnames = {'T2D','BMI','TG'};
npleio = 2;
plot_heatmap(est,tnames,npleio)
```

To generate heat maps of variance explained for each trait similar to the heat maps shown in our manuscript, one can use the function within the MODEL class, `h2plot`. Below is an example using the model output by the example script, `fit_model.m`:
```
traits = {'T2D','BMI','TG'};
cptnames = {'cpt1','cpt2'};
included_cpts = [4,5];
est.h2plot('traitNames',traits,'cptNames',cptnames,'whichCpts',included_cpts);
title('Variance explained by each component')
```

## Hypothesis testing
PDR has the ability to calculate p-values for testing a null model against an alternative model, where the null model contains a subset of the components in the alternative model (i.e., the null and alternative models are nested). The `test_cpts` function  within `MODEL` can be used to test the significance of each of the trait-specific components by iteratively removing these components one at a time and calculating a p-value. If the p-value is null, this indicates that the trait-specific component was needed to produce a better model fit. Otherwise, the p-value is significant. One can run the `test_cpts` function as follows:
```
pval = model.test_cpts(ecf);
```
Where `model` is the model that has been fit to the data, and `ecf` is the empirical characteristic function object that had been used to fit the model. See the documentation for more details on the other parameter options.

PDR can also give p-values for comparing models that the user specifies, as long as those models are nested. To do this, one would run the `test` function within `MODEL` as follows:
```
pval = nullmodel.test(altmodel,nullmodel_ecf);
```
Where `nullmodel_ecf` is the ECF that used to fit the null model.

Alternatively, one can calulate a p-value for a model against an alternative that represents the best that the model could fit the data. This gives a sense of how well the model fits the data, and whether there exists a larger model that could produce a better fit. This can be produced as follows:
```
pval = model.test(ecf);
```

## Calculating posterior mean effect sizes
Using the effect size distribution estimated by PDR, we can calculate posterior mean effect sizes as well as the membership scalars for each SNP-component pair. This can be done using the `predict` function within `MODEL` as follows:
```
[alpha, ~, posterior_scalars] = est.predict(data,'minWeight',1e-6);
```
Where `alpha` is the posterior mean effect sizes (# SNPs x # traits) and `posterior_scalars` are the posterior variances (# SNPs x # components) that indicate the probability that the SNP was drawn from a given component. `data` is the data object that was used for loading in the data (see "Loading Data"). `minWeight` is a parameter within `predict` which sets a threshold for discarding component combinations that are less likely. This helps with reducing the computational burden.

### Plotting posterior mean effect sizes
To plot scatter plots of the posterior mean effect sizes on pairs of traits similar to those in Fig. 2c-g and Fig. 3c-g, one can run the function `make_scatterplot` (located in `plotting`). Note that this function depends on `alpha` and `posterior_scalars`, which are outputs of the `predict` function (see above). You will also need to download `1kg_LD.HM3.window1cm.noblocks.mat`. This can be  done by running the following command:
```
wget https://www.dropbox.com/sh/mclm1urkxs8ga80/AABLDZRREAkGj5A1D3x8z_FOa/1kg_LD.HM3.window1cm.noblocks.mat?dl=0
```
Below is an example for plotting pairs of traits using the model output by `fit_model.m`:
```
cptnames = {'cpt1','cpt2'};
traits = {'T2D','BMI','TG'};
subplot_traits = [1 2; 1 3; 2 3];
included_cpts = [4,5];
load('/path/to/1kg_LD.HM3.window1cm.noblocks.mat','LDSNPs','RRb');
make_scatterplot(est,data,alpha,posterior_scalars,cptnames,traits,subplot_traits,included_cpts,LDSNPs,RRb)
```
This will make a scatter plot for every pair of traits specified in each row of `subplot_traits`.

## Calculating replication r<sup>2</sup>
To calculate expected and/or observed replication r<sup>2</sup>, one can run the function, `expected_observed_r2` (located in the subroutines directory). First, one must have the fitted PDR model (`est`) loaded, as well as the data (`data`) (see "Loading Data") and the posterior mean effect sizes (`alpha`) (see "Calculating posterior mean effect sizes"). Here, `traitidx` is 1 because this corresponds to T2D, the trait for which we are assessing replication r<sup>2</sup>. This function is run as follows:
```
traitnames = {'T2D','BMI','TG'};
replication_files = {'/path/to/replication_t2d_data.sumstats'};
traitidx = 1;
output = expected_observed_r2(est,data,traitnames,alpha,replication_files,traitidx);
```
Or, if there is no replication data and only predicted replication r<sup>2</sup> is desired, one can run:
```
output = expected_observed_r2(est,data,traitnames,alpha,{},traitidx);
```
Note that `output` is a struct with replication r<sup>2</sup> for PDR, MTAG, and the original summary statistics, and it will have either expected and observed, or only the expected values depending on whether replication_files is empty. `output` will also contain the numerical values for the per-SNP variance explained for both MTAG and PDR.

### Plotting replication r<sup>2</sup> results
A bar plot comparing the expected (and observed, if included) replication r<sup>2</sup> for the original summary statistics, MTAG, and PDR can be generated using the function, `r2_barplot` (located in the `plotting` directory).

The variance explained scatter plots comparing expected and observed variance can be generated using the `plot_varexp` function found within `plotting`. This plot can only be created if there is replication data available.

Both of these plotting methods require the output of `expected_observed_r2` (see above).

