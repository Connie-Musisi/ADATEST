#' Compute Test Statistic for Differential Abundance Analysis
#'
#' Computes a test statistic for each taxon based on the trained regression parameters (scores) from the PLS model.
#' Samples are grouped, sorted by abundance, and aligned based on group-level means to reduce noise and enhance interpretability.
#' This function is used internally by \code{Orig_Alz_Perm()}.
#'
#' @param org_otu A numeric matrix representing the scaled OTU table from the original dataset, with samples in rows and taxa in columns.
#' @param Train_parest A numeric vector of trained regression parameters (scores) obtained from \code{Train_analyze()}.
#' @param n1 Integer. Number of samples in group 1. Passed from \code{Orig_Alz_Perm()}.
#' @param n2 Integer. Number of samples in group 2. Passed from \code{Orig_Alz_Perm()}.
#' @param B Integer. Number of permutations (included for interface consistency; not used directly in this function).
#'
#' @return A numeric vector of test statistics, one for each taxon.
#'
#' @details
#' For each taxon:
#' \itemize{
#'   \item The function compares group-specific mean abundances to determine ordering.
#'   \item OTUs are sorted within each group to reduce intra-group variability.
#'   \item Taxa are then aligned by group and multiplied by the regression scores to yield a final test statistic.
#' }
#' The output is used to assess significance against permuted null distributions in the overall ADATEST pipeline.
#'
#' @keywords internal
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
