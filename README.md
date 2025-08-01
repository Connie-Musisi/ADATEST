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
All data supporting the findings of this study are publicly available in this repository:

- **Simulated Data**  
  Available under [`SimData/`](./SimData):
  - [`SimData/IBD/`](./SimData/IBD): Simulations for the IBD setting under multiple scenarios.
  - [`SimData/Dietswap/`](./SimData/Dietswap): Simulations for the DietSwap setting.
  - [`SimData/simulations R code/`](./SimData/simulations%20R%20code): R scripts used to generate the simulated datasets.

- **Analysis Code**  
  See the [`AnalysisCode/`](./AnalysisCode) folder for all scripts used in pre-processing, model training, and test statistic evaluation.

- **Results**  
  The [`Results/`](./Results) folder contains Excel files with processed outputs, including test statistics and FDR-adjusted p-values.

## Example
See example vignette for more explanation on how to use ADATEST on a real dataset.

## License
This project is licensed under the MIT License. See the [`LICENSE/`](ADATEST/LICENSE) file for details.



This repository ensures reproducibility of all analyses presented in the manuscript. No restrictions apply to data access.

