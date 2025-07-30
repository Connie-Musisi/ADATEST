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
## Code
See example vignette for more explanation on how to use ADATEST on a real dataset.

## Data Aavilability
The simulated datasets can be found in the  folder SimData
