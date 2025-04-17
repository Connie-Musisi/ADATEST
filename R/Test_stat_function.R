#' Compute Test Statistic for Differential Abundance Analysis
#'
#' This function calculates a test statistic for each taxon using a PLS estimates/scores from the trained model.
#' It sorts samples within each group and aligns them according to their mean abundance before computing the test statistic.
#'
#' @param org_otu A matrix representing the OTU table with samples in rows (OTU tavle from original dataset)
#' @param Train_parest A vector of estimated parameters (scores) used for computing the test statistic.
#' @param n1 The number of samples in group 1 (default: half of the dataset).
#' @param n2 The number of samples in group 2 (default: half of the dataset).
#' @param B Number of permutations used for significance testing (default: 100).
#'
#' @return A vector of computed test statistics for each taxon.
#'
#' @export
Test_Statistic <- function(org_otu, Train_parest, n1, n2, B) {
  n.taxa <- ncol(org_otu)
  
  otu1 <- org_otu[1:n1, ]
  otu2 <- org_otu[(n1 + 1):(n1 + n2), ]
  
  taxa.means1 <- colMeans(otu1)
  taxa.means2 <- colMeans(otu2)
  
  # Identify group with smaller mean
  smallest <- ifelse(taxa.means1 - taxa.means2 < 0, 1, 2) 
  
  # Sort within groups
  otu1.sorted <- t(apply(otu1, 2, sort, decreasing = FALSE))
  otu2.sorted <- t(apply(otu2, 2, sort, decreasing = FALSE))
  
  # Align groups based on mean abundance
  select1 <- (smallest == 1)  
  select2 <- (smallest == 2)  
  otu.final <- matrix(nrow = n.taxa, ncol = (n1 + n2))
  otu.final[select1, (n1 + 1):(n1 + n2)] <- otu2.sorted[select1, ]
  otu.final[select1, 1:n1] <- otu1.sorted[select1, ]
  otu.final[select2, (n1 + 1):(n1 + n2)] <- otu1.sorted[select2, 1:n2]
  otu.final[select2, 1:n2] <- otu2.sorted[select2, ]
  
  # Compute test statistic using estimated parameters
  org_pred <- otu.final %*% Train_parest
  
  return(stats = org_pred)
}
