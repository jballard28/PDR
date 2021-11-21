# PDR

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

There are many parameters in this script that can be altered according to which type of model you would like to fit. For example, the number of rank-one ('factor-like') and full-rank ('generic pleiotropic') components can be altered. If more traits or more components are being fit, we recommend increasing the number of initializations (`ninit`) and also the maximum number of gradient descent steps (`gdsteps`). If the code is giving warnings that there is an insufficient number of sampling times compared to the number of parameters, we recommend reducing the `rotation_tolerance` parameter in the `initialize_fit` function.

## Visualizing the components as heat maps
See `example_scripts/plot_heatmaps.m` for an example using the output of `fit_model.m`.

This script visualizes the component covariance matrices normalized such that the total variance explained across all components for each trait is equal to 1.

To generate heat maps of variance explained for each trait similar to the heat maps shown in our manuscript, see `plot_varexp_heatmaps.m`. To generate the latter heat maps, you will need to download `1kg_LD.HM3.window1cm.noblocks.mat`. This can be  done by running the command:
```
wget https://www.dropbox.com/sh/mclm1urkxs8ga80/AABLDZRREAkGj5A1D3x8z_FOa/1kg_LD.HM3.window1cm.noblocks.mat?dl=0
```

## Hypothesis testing

## Calculating posterior mean effect sizes

## Calculating replication r<sup>2</sup>
