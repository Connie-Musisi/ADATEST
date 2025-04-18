#' Adaptive Test for Differential Abundance Analysis 
#'
#' This function performs an adaptive differential abundance test (ADATEST) using a permutation-based approach.
#' It applies total sum scaling (TSS), trims the dataset to remove taxa that are not prevalent,
#' subsets the data to create the Pseudo dataset and finally the training dataset,
#' performs parameter estimation, conducts a permutation test, and evaluates performance metrics.
#'
#' @param physeq A phyloseq object containing microbiome OTU/ASV count data.
#' @param group_var Name of the grouping variable in sample metadata (default: "group").
#' @param group_levels A character vector of two levels in \code{group_var} to be compared (default: c("0", "1")).
#' @param t0 Lower threshold for effect size classification (default: 1).
#' @param t1 Upper threshold for moderate effect classification (default: 1.5).
#' @param t2 Maximum bound for strong effect classification (default: +Inf).
#' @param B Number of permutations for significance testing (default: 1000).
#' @param empirical_adjust Logical. If TRUE, applies empirical adjustment for unequal group sizes (default: TRUE).
#' @param p.adjust.method Method for p-value adjustment; options include "BH", "bonferroni", "fdr", etc. (default: "BH").
#' @param n.taxa0 Number of non-differentially abundant taxa, if known. Used for tuning (default: NULL).
#' @param seed Random seed for reproducibility (default: \code{set.seed(19)}).
#'
#' @return A list containing:
#' \item{Original_results}{Results from the permutation test, including p-values, adjusted p-values, test statistics, and DA status.}
#' \item{Unscaled_data}{A list with the original phyloseq object and the training dataset used.}
#' \item{Train_results}{Results from the training phase including parameter estimates and prediction performance.}
#' \item{Pseudo}{Pseudo data summaries, including delta_med and centered log fold changes.}
#'
#' @export
ADATEST <- function(physeq, group_var = "group",
                    group_levels = c("0", "1"),
                    t0=1, t1=1.5, t2=+Inf, B=1000,
                    empirical_adjust = TRUE, 
                    p.adjust.method = "BH",
                    n.taxa0=NULL,
                    seed=set.seed(19)){
  # Total Sum Scaling
  simdata_filter <- transform_sample_counts(physeq, function(x) { x / sum(x)})
  
  # Trimming
  simdata_filter <- signtrans::Trim(simdata_filter,minReads = 0.001, minPrev = 0.20) 
  
  simdata_filter <- transform_sample_counts(simdata_filter, function(x) { x / sum(x)})  
 
  cat("running DATA")
  # Subset data for training
  Pseudo_Train <- PseudoData(simdata_filter, group_var = group_var, 
                             group_levels = group_levels,
                             t0 = t0, t1, n.taxa0 = n.taxa0, t2, 
                             empirical_adjust = empirical_adjust)
  
  Train <- Pseudo_Train$Train
  
  Train_nan <- NA.remove(Train)
  Train_res <- Train_analyze(Train_nan, seed=seed) 
  
  # Permutation Test
  cat("running permutation test")
  Orig_results <- Orig_Alz_Perm(simdata_filter, group_var = group_var,
                                group_levels = group_levels, 
                                Train_parest = Train_res$Train_parest, 
                                p.adjust.method = p.adjust.method, B = B)
  
  
  return(list(Original_results = Orig_results,
              Unscaled_data = list(Train_data = Train_nan, Original_data = simdata_filter),
              Train_results = Train_res,
              Pseudo = list(Pseudo_Train$delta_med, Pseudo_Train$lfc_centered, Pseudo_Train$pseudo_lfc)
              ))
}