#' Adaptive Test for Differential Abundance Analysis 
#'
#' This function performs an adaptive differential abundance test using a permutation-based approach.
#' It applies total sum scaling (TSS), trims the dataset to remove taxa that are not prevalent,
#' subsets the data to create the Pseudo dataset and finally the training dataset,
#' performs parameter estimation, conducts a permutation test, and evaluates performance metrics.
#'
#' @param physeq A phyloseq object containing OTU data.
#' @param t0 Lower bound for data transformation (default: 1).
#' @param t1 Upper bound for data transformation (default: 1.5).
#' @param t2 Maximum bound for data transformation (default: +Inf).
#' @param B Number of permutations used for significance testing (default: 1000).
#'
#' @return A list containing:
#'   - `Original_results`: Results from the permutation test.
#'   - `Unscaled_data`: List containing training data and original data.
#'   - `Train_results`: Results from training data analysis.
#'   - `Results`: Evaluation metrics including FDR, sensitivity, specificity, and power.
#'
#' @export
ADATEST <- function(physeq, t0=1, t1=1.5, t2=+Inf, B=1000){
  # Total Sum Scaling
  simdata_filter <- transform_sample_counts(physeq, function(x) { x / sum(x)})
  
  # Trimming
  simdata_filter <- signtrans::Trim(simdata_filter,minReads = 0.001, minPrev = 0.20) 
  
  simdata_filter <- transform_sample_counts(simdata_filter, function(x) { x / sum(x)})  
  tmp <- tax_table(simdata_filter)
  n.taxa <- nrow(tmp)
  n.taxa1 <- sum(as.data.frame(tmp)$isDA == "TRUE")
  
  cat("running DATA")
  # Subset data for training
  Pseudo_Train <- PseudoData(simdata_filter, t0, t1, n.taxa0 = n.taxa - n.taxa1, t2, empirical_adjust = FALSE)
  
  Train <- Pseudo_Train$Train
  
  Train_nan <- NA.remove(Train)
  Train_res <- Train_analyze(Train_nan) 
  
  # Permutation Test
  cat("running permutation test")
  Orig_results <- Orig_Alz_Perm(simdata_filter, Train_parest = Train_res$Train_parest, 
                                p.adjust.method = "qvalue", B = B)
  
  # Evaluation
  Results <- eval(Orig_results$Original_results)
  
  return(list(Original_results = Orig_results,
              Unscaled_data = list(Train_data = Train_nan, Original_data = simdata_filter),
              Train_results = Train_res,
              Results = Results))
}