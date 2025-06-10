#' Adaptive Test for Differential Abundance Analysis 
#'
#' This function performs an adaptive differential abundance test (ADATEST) using a permutation-based approach.
#' It applies total sum scaling (TSS), trims the dataset to remove taxa that are not prevalent, reaarranges the groups
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
#' @return A named \code{list} with components:
#' \describe{
#'   \item{Original_results}{\code{data.frame} of the permutation test output on the original dataset,
#'     containing columns:
#'     \describe{
#'       \item{Taxon}{feature identifier}
#'       \item{Stat_Obs}{observed test statistic}
#'       \item{p}{raw permutation p-value}
#'       \item{p.adjusted}{adjusted p-value via \code{p.adjust.method}}
#'       \item{DA}{logical indicator of differential abundance at the chosen FDR level}
#'     }
#'   }
#'   \item{Train_results}{\code{list} containing:
#'     \describe{
#'       \item{Train_parest}{estimated score parameters from the training phase (PLS coefficients)}
#'       \item{Train_performance}{metrics (e.g., accuracy, AUC) evaluating predictive fit on pseudo-data}
#'       \item{Train_scaled}{rescaled training matrix used for scoring}}
#'   }
#'   \item{Unscaled_data}{\code{list} with:
#'     \describe{
#'       \item{Original_data}{the trimmed-and-TSS-normalized \code{phyloseq} object passed to training}
#'       \item{Train_data}{pseudo-taxa training data before scaling or permutation}}
#'   }
#'   \item{Scaled_data}{\code{list} with:
#'     \describe{
#'       \item{Original_scaled}{rescaled counts of the original dataset used in test statistic computation}
#'       \item{Train_scaled}{rescaled counts of pseudo-data used for training}}
#'   }
#'   \item{Pseudo_results}{output from \code{PseudoData()}, including pseudo-taxa summaries and thresholds}
#'   \item{p_value}{numeric vector of raw p-values (extracted from \code{Original_results}$p)}
#'   \item{p_adj}{numeric vector of adjusted p-values (extracted from \code{Original_results}$p.adjusted)}
#'   \item{test_statistic}{numeric vector of observed test statistics (extracted from \code{Original_results}$Stat_Obs)}
#' }
#'  \item  {data_Adjust2Median}{Data corrected for compositionality} 
#'
#' @details
#' ADATEST adapts to dataset-specific characteristics via a three-step procedure:
#' 1. **Trimming & Scaling**: Applies TSS and filters low-prevalence taxa to mitigate sparsity.
#' 2. **Pseudo-Data Generation & Training**: Constructs pseudo-taxa from the larger group, assigns labels
#'    based on effect-size thresholds (\code{t0}, \code{t1}, \code{t2}), and fits a PLS model to estimate
#'    adaptive scores.
#' 3. **Permutation Testing**: Computes test statistics on the real data using estimated scores,
#'    permutes sample labels \code{B} times, then derives p-values and controls FDR via \code{p.adjust}. 
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
  simdata_tss <- transform_sample_counts(physeq, function(x) { x / sum(x)})
  
  # Trimming
  simdata_filter <- signtrans::Trim(simdata_tss,minReads = 0.001, minPrev = 0.20) 
  
  # Renormalize after trimming
  simdata_tss2 <- transform_sample_counts(simdata_filter, function(x) { x / sum(x)})  

   # Adjust to median (your custom function)
  simdata_adjusted <- Adjust2Median(simdata_tss2)
  
  print(simdata_adjusted)

  
  # Reorder data so that the group with more samples comes first (n1>n2)
  samdata <- as.data.frame(sample_data(simdata_adjusted))
  grp_vec <- samdata[[group_var]]
  counts <- table(grp_vec)
  if (counts[group_levels[1]] >= counts[group_levels[2]]) {
    primary <- group_levels[1]
    secondary <- group_levels[2]
  } else {
    primary <- group_levels[2]
    secondary <- group_levels[1]
  }
  
  ordered_samps <- c(rownames(samdata)[grp_vec == primary],
                     rownames(samdata)[grp_vec == secondary])
  samdata <- samdata[ordered_samps, , drop = FALSE]
  sample_data(simdata_adjusted) <- sample_data(samdata)
  
  if (taxa_are_rows(simdata_adjusted)) {
    otu_mat <- as(otu_table(simdata_adjusted), "matrix")
    otu_mat <- otu_mat[, ordered_samps]
    otu_table(simdata_adjusted) <- otu_table(otu_mat, taxa_are_rows = TRUE)
  } else {
    otu_mat <- as(otu_table(simdata_adjusted), "matrix")
    otu_mat <- otu_mat[ordered_samps, ]
    otu_table(simdata_adjusted) <- otu_table(otu_mat, taxa_are_rows = FALSE)
  }
  
  cat("running DATA")
  # Subset data for training
  Pseudo_Train <- PseudoData(simdata_adjusted, group_var = group_var, 
                             group_levels = c(primary,secondary),
                             t0 = t0, t1, n.taxa0 = n.taxa0, t2, 
                             empirical_adjust = empirical_adjust)
  
  Train <- Pseudo_Train$Train
  
  Train_nan <- NA.remove(Train)
  Train_res <- TrainAlz(Train_nan,group_var = "new_group", 
                             group_levels = c("0","1"), seed=seed) 
  
  # Permutation Test
  cat("running permutation test")
  Orig_results <- OrigAlz(simdata_adjusted, group_var = group_var,
                                group_levels = c(primary,secondary), 
                                Train_parest = Train_res$Train_parest, 
                                p.adjust.method = p.adjust.method, B = B)
  
  return(list(Original_results = Orig_results,
              Train_results = Train_res,
              Pseudo_results = Pseudo_Train,
              Unscaled_data = list(Train_data = Train_nan, Original_data = simdata_tss2),
              Scaled_data = list(Train_scaled=Train_res$Train_scaled, Original_scaled=Orig_results$org_scaled),
              data_Adjust2Median = simdata_adjusted, 
              p_value = Orig_results$p,
              p_adj = Orig_results$p.adjusted,
              test_statistic = Orig_results$Stat_Obs
              
  ))
}
