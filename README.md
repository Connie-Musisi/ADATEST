# ADATEST: An Adaptive Test for Differential Abundance in Microbiome Studies
## Overview
The ADATEST package implements a novel adaptive test for differential abundance analysis (DAA) using order statistics in microbiome studies when there are two groups. This method is designed to handle the complex characteristics of microbiome data, such as sparsity, compositionality, and overdispersion.

### Key Steps
The ADATEST algorithm involves the following three steps:
1. **Correction for compositionality:** The relative abundances of the original dataset are centered using the median log-fold change to remove compositional bias.
2. **Construct a pseudo-dataset:** Generate all unique taxa pairs from one group and use these to build a labeled training dataset for model fitting.
3. **Estimate scores from training data:** Fit a linear model (Partial Least Squares) to the training dataset to learn taxon-level importance scores.
4. **Compute test statistics on the original data:** Apply the estimated scores to the original dataset to compute test statistics and the p-values obtained from the permutation null distribution of the test statistics are adjusted to control the FDR.

## Installation
To install the latest version of ADATEST from GitHub, you can use the 'devtools' package.
```r
install.packages("devtools")
library(devtools)
devtools::install_github("Connie-Musisi/ADATEST")
```

### Required packages
The packages required to run ADATEST are `phyloseq`, `signtrans`, `qvalue` and `stats`.
    
## Data Availability
All data supporting the findings of this study are publicly available in this repository. The simulated datasets are available under [`SimData/`](./SimData): \
  `IBD data` and `Dietswap data` are the simulated datasets used in analysis for the results presented in the main manuscript. For each source data, there are 12 different settings. \
  `IBD_NoComp data`, `Dietswap_NoComp data` were used for the no compositionality analysis which is presented in the Supplementary file. The results from `NegativeBinomial data` and `SPsimSeq data` simulated datasets are aslo presneted in the Supplementary. \
  `DATA_SIMULATION_FRAMEWORK.pptx` gives an overview of the simulated framework settings. 

## Analysis Code 
  See the [`All_analysis_Code/`](./All_analysis_Code) for all the R scripts used in analysis. \
  The folder [`Code_for_analysis_results_manuscript/`](./All_analysis_Code/Code_for_analysis_results_manuscript) contains:
  1. The file `Functions_for_methods` is where a function is made for each method to output the test statistics, raw and adjusted p-values. It also contains the `eval` function which calculates the FDR, Sensitivity and Type I error.
  2. The file `Analysis_for_all_methods` uses all the functions of the different methods to do analysis on each data setting. The output is a dataframe containing all the FDR and Sensitivity values for all the methods for all the datasets. 
  3. The file `Code_for_plotting_graphs` contains the R code for all the plots that are presented in the main manuscript and Supplementary file.
 The folder [`Code_for_simulating_data/`](./All_analysis_Code/Code_for_simulating_data) provides the R scripts `Simulate_data` used to generate the simulated datasets, starting with using `MIDASim` to simulate the counts and then creating a phyloseq object for each simulated count dataset.

  
## Results  
  The [`Results/`](./Results) folder contains Excel files with processed outputs, including FDR and Sensitivity values for all the methods. \
  It also contains the .RData reults from ADATEST analysis, this includes the unscaled and scaled original and training datasets, score values, test statistics, raw and adjusted p-values. \
  The output is:
  | Value            | Description                                   |
|------------------|-----------------------------------------------|
| `Original_results` | Results from the `OrigAlz` including the scaled Original data, test statistics, raw and adjusted p-values and the new differential abundance indicator |
| `Train_results` | Results from the `TrainAlz` including the scaled Training dataset used to fit the model, predictions on the training data and the parameter estimates from the model which are the sccores |
| `Pseudo_results` | Results from the `PseudoData` funcrtion including the Training dataset used as input in `TrainAlz` and   |
| `Unscaled_data`    | Original and Training datasets before the `Normalize` function is applied   |
| `Scaled_data`   | Scaled Original and training datasets using the `Normalize` function       |
| `data_Adjust2Median` | The dataset adjusted to correct for compositionality usind the `Adjust2Median` function    |
| `p_value` | Raw p-values |
| `p_adj` | Adjusted p-values   |
| `test_statistic` | Test statistic values for each taxon |


## Example
The input for `ADATEST` is a phyloseq object with an OTU table, meta data and a taxa table. See example [`vignette`](./vignettes/Example.Rmd) for more explanation on how to use `ADATEST` on a real dataset.

## License
This project is licensed under the MIT License. See the [`LICENSE/`](./LICENSE) file for details. \
This repository ensures reproducibility of all analyses presented in the manuscript. No restrictions apply to data access.

## References
[1] He, M., Zhao, N. & Satten, G.A. 
MIDASim: a fast and simple simulator for realistic microbiome data. *Microbiome* 12, 135 (2024). 
[https://doi.org/10.1186/s40168-024-01822-z](https://doi.org/10.1186/s40168-024-01822-z)

[2] O’Keefe, S., Li, J., Lahti, L. et al. 
Fat, fibre and cancer risk in African Americans and rural Africans. *Nat Commun* 6, 6342 (2015). 
[https://doi.org/10.1038/ncomms7342](https://doi.org/10.1038/ncomms7342)

[3] Lloyd-Price, J., Arze, C., Ananthakrishnan, A.N. *et al.* 
Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. *Nature*, 569, 655–662 (2019). 
[https://doi.org/10.1038/s41586-019-1237-9](https://doi.org/10.1038/s41586-019-1237-9)

[4]  Kodalci L, Thas O (2023) 
Simple and flexible sign and rank-based methods for testing for differential abundance in microbiome studies. *PLoS ONE* 18(9): e0292055. 
[https://doi.org/10.1371/journal.pone.0292055](https://doi.org/10.1371/journal.pone.0292055)
