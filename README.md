# ADATEST: An Adaptive Test for Differential Abundance in Microbiome Studies
## Overview
The ADATEST package implements a novel adaptive test for differential abundance analysis (DAA) using order statistics in microbiome studies when there are two groups. This method is designed to handle the complex characteristics of microbiome data, such as sparsity, compositionality, and overdispersion. It involves three key steps: (1) constructing a pseudo-dataset by generating all unique taxa pairs from one group, and then obtaining a training dataset, (2) estimating scores from the training dataset using a linear model (Partial Least Squares), and (3) using the score estimates to compute the test statistics of the original dataset to detect differential abundance.

## Installation
To install the latest version of ADATEST from GitHub, you can use the 'devtools' package.
```r
install.packages("devtools")
library(devtools)
devtools::install_github("Connie-Musisi/ADATEST")
```

## Data Aavilability
All data supporting the findings of this study are publicly available in this repository. The simulated datasets are available under [`SimData/`](./SimData): \
  `IBD data` and `Dietswap data` are the simulated datasets used in analysis for the results presented in the main manuscript. For each source data, there are 12 different settings. \
  `IBD_NoComp data`, `Dietswap_NoComp data` were used for the no compositionality analysis which is presented in the Supplementary file. The results from `NegativeBinomial data` and `SPsimSeq data` simulated datasets are aslo presneted in the Supplementary. \
   R scripts used to generate the simulated datasets.

## Analysis Code 
  See the [`Analysis Code/`](./AnalysisCode) folder for all the R scripts used evaluation of ADATEST and all competitor methods.
  The file `Functions_for_methods` is where a function is made for each method to output the test statistics, raw and adjusted p-values. It also contains the `eval` function which calculates the FDR, Sensitivity and Type I error. \
  The file `Analysis_for_methods` uses all the functions of the different methods to do analysis on each data setting. The output is a dataframe containing all the FDR and Sensitivity values for all the methods for all the datasets. \
  The file `Code_for_plotting_graphs` contains the R code for all the plots that are presented in the main manuscript and Supplementary file.
  
## Results  
  The [`Results/`](./Results) folder contains Excel files with processed outputs, including FDR and Sensitivity values for all the methods. \
  It also contains the .RData reults from ADATEST analysis, this includes the unscaled and scaled original and training datasets, score values, test statistics, raw and adjusted p-values. \
  The output is
  | Value            | Description                                   |
|------------------|-----------------------------------------------|
| Original_results | tetrachoric correlation (presence/absence)    |
| Train_results | Pearson correlation among non-zeros         |
| Pseudo_results | proportion of non-zeros per taxon             |
| Unscaled_data    | proportion of non-zeros per sample            |
| Scaled_data   | mean relative abundance per taxon             |
| data_Adjust2Median      | vector of non-zero abundances (each taxon)    |
| p_value | generalized gamma parameters (parametric mode only) |
| p_adj |  |
| test_statistic |  |


## Example
See example [`vignette`](./vignettes/Example.Rmd) for more explanation on how to use ADATEST on a real dataset.

## License
This project is licensed under the MIT License. See the [`LICENSE/`](./LICENSE) file for details.



This repository ensures reproducibility of all analyses presented in the manuscript. No restrictions apply to data access.

