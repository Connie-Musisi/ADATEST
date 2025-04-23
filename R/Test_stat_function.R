#' Compute Test Statistic for Differential Abundance Analysis
#'
#' Computes a test statistic for each taxon based on the trained regression parameters (scores) from the PLS model.
#' Samples are grouped, sorted by abundance, and aligned based on group-level means to reduce noise and enhance interpretability.
#' This function is used internally by \code{Orig_Alz_Perm()}.
#'
#' @param otu A numeric matrix representing the scaled OTU table from the original dataset, with samples in rows and taxa in columns.
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

Test_Statistic <- function(otu, Train_parest, n1, n2, B) {
  n.taxa <- ncol(otu)
  otu1 <- otu[1:n1, , drop = FALSE]
  otu2 <- otu[(n1 + 1):(n1 + n2), , drop = FALSE]
  
  taxa.means1 <- colMeans(otu1)
  taxa.means2 <- colMeans(otu2)
  smallest <- ifelse(taxa.means1 - taxa.means2 < 0, 1, 2)
  
  # Pad the sorted data to equal lengths (max(n1, n2)) using last value
  pad_sort <- function(x, len) {
    x_sorted <- sort(x)
    if (length(x_sorted) < len) {
      pad <- rep(tail(x_sorted, 1), len - length(x_sorted))
      return(c(x_sorted, pad))
    } else {
      return(x_sorted)
    }
  }
  
  max_n <- max(n1, n2)
  otu1.sorted <- t(apply(otu1, 2, pad_sort, len = max_n))
  otu2.sorted <- t(apply(otu2, 2, pad_sort, len = max_n))
  
  # Now fill final matrix of size: n.taxa × (n1 + n2)
  otu.final <- matrix(nrow = n.taxa, ncol = (n1 + n2))
  
  # For taxa where group 1 is smaller
  sel1 <- smallest == 1
  otu.final[sel1, 1:n1] <- otu1.sorted[sel1, 1:n1, drop = FALSE]
  otu.final[sel1, (n1 + 1):(n1 + n2)] <- otu2.sorted[sel1, 1:n2, drop = FALSE]
  
  # For taxa where group 2 is smaller
  sel2 <- smallest == 2
  otu.final[sel2, 1:n2] <- otu2.sorted[sel2, 1:n2, drop = FALSE]
  otu.final[sel2, (n2 + 1):(n2 + n1)] <- otu1.sorted[sel2, 1:n1, drop = FALSE]
  
  # Apply test statistic
  org_pred <- otu.final %*% Train_parest
  
  return(org_pred)
}
