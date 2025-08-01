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

Simulated data: The simulated Data is located in the SimData/IBD and SimData/Dietswap folder for the simulation outputs for both the IBD and DietSwap settings respectively, under multiple scenarios. The R code used for simulating the data can be found in SimData/simulations R code

Analysis: Full R scripts used to generate results can be found in the AnalysisCode/ folder.

Results: Results: Processed results, including test statistics and FDR-adjusted p-values, are provided in Excel format in the Results/ directory.

## Example
See example vignette for more explanation on how to use ADATEST on a real dataset.

## License
This project is licensed under the MIT License. See the LICENSE file for details.



This repository ensures reproducibility of all analyses presented in the manuscript. No restrictions apply to data access.

