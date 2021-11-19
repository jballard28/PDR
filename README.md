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

Download the summary statistics that will be input to PDR. For some example data, one can run the commands:\
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
weightsfp = '/path/to/1000G_Phase3_weights_hm3_no_MHC_unzipped/weights.hm3.noMHC.';

data = DATA(clustertraits, 'LDscoreDir', l2fp, 'WeightsDir', weightsfp);
```
