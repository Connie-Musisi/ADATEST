#' Construct Pseudo Training Dataset for DA Testing (Internal)
#'
#' Constructs a pseudo-dataset from microbiome data to support adaptive differential abundance (DA) testing.
#' This internal function is called by \code{ADATEST()} to generate synthetic pseudo-taxa for model training, 
#' based on taxon pairing, log fold change (LFC) classification, and optional empirical adjustment.
#'
#' @param data A \code{phyloseq} object containing filtered microbiome data.
#' @param group_var Character. The name of the grouping variable in the sample metadata. Passed from \code{ADATEST()}.
#' @param group_levels Character vector of length two. The two levels of \code{group_var} to be compared. Passed from \code{ADATEST()}.
#' @param t0 Numeric. Lower threshold for classifying non-differential pseudo-taxa. Passed from \code{ADATEST()}.
#' @param t1 Numeric. Threshold for identifying moderate differential abundance. Passed from \code{ADATEST()}.
#' @param t2 Numeric. Upper bound for classifying pseudo-taxa. Passed from \code{ADATEST()}.
#' @param empirical_adjust Logical. If \code{TRUE}, downsamples group 2 using an empirical distribution. Passed from \code{ADATEST()}.
#' @param n.taxa0 Optional. Number of non-DA pseudo-taxa to retain. Passed from \code{ADATEST()}.
#'
#' @return A list with components:
#' \describe{
#'   \item{pseudo_lfc}{Vector of log fold changes for pseudo-taxa.}
#'   \item{Train}{A \code{phyloseq} object containing the training dataset with pseudo-taxa labeled as \code{I_0} or \code{I_1}.}
#' }
#'
#' @details
#' This function is not intended to be called directly by users. It is part of the internal ADATEST pipeline:
#' \enumerate{
#'   \item Taxa are paired within the reference group to generate pseudo-taxa.
#'   \item Pseudo-taxa are filtered for prevalence and abundance.
#'   \item LFCs are computed and centered; taxa are classified into \code{I_0}, \code{I_1}, or \code{I_B}.
#'   \item The training set includes only pseudo-taxa from \code{I_0} and \code{I_1}.
#'   \item If needed, empirical resampling balances group sizes before model training.
#' }
#'
#' @keywords internal
#' @import phyloseq
#' @export

PseudoData <- function(simdata_filter,
                       group_var,
                       group_levels,
                       t0, 
                       t1, 
                       t2,
                       empirical_adjust,
                       n.taxa0) {
  
  samdata <- as.data.frame(sample_data(simdata_filter))
  grp_vec <- samdata[[group_var]]
  
  # Subset samples for each group
  samps_grp1 <- rownames(samdata)[grp_vec == group_levels[1]]
  samps_grp2 <- rownames(samdata)[grp_vec == group_levels[2]]
  
  simdata_group1 <- prune_samples(samps_grp1, simdata_filter)
  simdata_group2 <- prune_samples(samps_grp2, simdata_filter)
  
  otutab_group1 <- if (taxa_are_rows(simdata_group1)) t(otu_table(simdata_group1)) else otu_table(simdata_group1)
  otutab_group2 <- if (taxa_are_rows(simdata_group2)) t(otu_table(simdata_group2)) else otu_table(simdata_group2)
  
  n_samples_group1 <- nrow(otutab_group1)
  n_samples_group2 <- nrow(otutab_group2)
  n.taxa <- ncol(otutab_group1)
  if (is.null(n.taxa0)) {
    n.taxa0 <- round(n.taxa / 2)
  }
  
  taxon_names <- taxa_names(simdata_group1)
  taxon_name_pairs <- c()
  
  for (i in 1:(n.taxa - 1)) {
    use <- nrow(otutab_group1)
    pairotu_subgroup1 <- t(matrix(otutab_group1[, i], nrow = n.taxa - i, ncol = use, byrow = TRUE))
    pairotu_subgroup2 <- otutab_group1[, (i + 1):n.taxa]
    
    taxon_name_pairs <- c(taxon_name_pairs, paste(rep(taxon_names[i], n.taxa - i), taxon_names[(i + 1):n.taxa], sep = "-"))
    
    if (i == 1) {
      combined_paired_otu <- rbind(pairotu_subgroup1, pairotu_subgroup2)
    } else {
      combined_paired_otu <- cbind(combined_paired_otu, rbind(pairotu_subgroup1, pairotu_subgroup2))
    }
  }
  
  # Compute log-fold change for each pseudo-taxon
  lfc_pairs <- apply(combined_paired_otu, 2, function(col) {
    g1 <- col[1:n_samples_group1]
    g2 <- col[(n_samples_group1 + 1):(n_samples_group1 + n_samples_group2)]
    min_nz <- min(col[col != 0], na.rm = TRUE)
    log2((mean(g2) + min_nz) / (mean(g1) + min_nz))
  })
  
  # Keep all pseudo-taxa (or optionally filter by abs(LFC) > threshold)
  idx.keep.freq0 <- apply(combined_paired_otu, 2, function(x) mean(x == 0)) < 0.8
  idx.keep.avg   <- apply(combined_paired_otu, 2, function(x) mean(x)) > 0.001
  idx.keep <- idx.keep.freq0 & idx.keep.avg
  
  combined_paired_otu <- combined_paired_otu[, idx.keep]
  lfc_pairs <- lfc_pairs[idx.keep]
  taxon_name_pairs <- taxon_name_pairs[idx.keep]
  
  # Construct phyloseq object for training
  tax_data <- tax_table(as.matrix(taxon_name_pairs))
  rownames(tax_data) <- taxon_name_pairs
  colnames(tax_data) <- "Pseudo_taxon"
  
  samdata_new <- rbind(sample_data(simdata_group1), sample_data(simdata_group2))
  samdata_new <- sample_data(as.data.frame(samdata_new))
  
  otu_table_new <- otu_table(combined_paired_otu, taxa_are_rows = FALSE)
  colnames(otu_table_new) <- taxon_name_pairs
  rownames(otu_table_new) <- rownames(samdata_new)
  
  training_data <- phyloseq(otu_table_new, tax_data, samdata_new)
  
  # Optional: downsample group 2 via empirical adjustment
  if (empirical_adjust) {
    train_group2 <- otu_table_new[(n_samples_group1 + 1):(2 * n_samples_group1), ]
    redtrain_group2 <- apply(train_group2, 2, function(x) {
      set.seed(19)
      emp_dist <- ecdf(x)
      quantile(emp_dist, probs = runif(n_samples_group2))
    })
    train_otu_table <- otu_table(rbind(otu_table_new[1:n_samples_group1, ],
                                       redtrain_group2), taxa_are_rows = FALSE)
    rownames(train_otu_table) <- rownames(simdata_filter@sam_data)
    train_sample_data <- sample_data(simdata_filter)[, group_var, drop = FALSE]
    training_data <- phyloseq(train_otu_table, tax_data, train_sample_data)
  }
  
  # Compute log-fold change again (on final training data)
  otu_train_filter <- training_data@otu_table
  data_group1 <- otu_train_filter[1:n_samples_group1, ]
  data_group2 <- otu_train_filter[(n_samples_group1 + 1):(n_samples_group1 + n_samples_group2), ]
  
  mean_group1 <- colMeans(data_group1)
  mean_group2 <- colMeans(data_group2)
  min.without.zero <- apply(otu_train_filter, 2, function(x) min(x[x != 0]))
  log_fold_change <- log2((mean_group2 + min.without.zero) / (mean_group1 + min.without.zero))
  
  # Label pseudo-taxa based on raw LFC
  tax_df <- as.data.frame(tax_table(training_data))
  tax_df$lfc <- log_fold_change
  tax_df$group_ind <- ifelse(abs(log_fold_change) < t0, "I_0",
                             ifelse(abs(log_fold_change) > t1 & abs(log_fold_change) < t2, "I_1", "I_B"))
  
  # Add group info to sample_data
  sample_df <- as.data.frame(sample_data(training_data))
  sample_df$new_group <- c(rep(0, n_samples_group1), rep(1, n_samples_group2))
  sample_data_new <- sample_data(sample_df)
  
  train_final <- phyloseq(otu_table(training_data), sample_data_new, tax_table(as.matrix(tax_df)))
  
  # Keep only I0 and I1 pseudo-taxa
  taxa_I0 <- tax_df$Pseudo_taxon[tax_df$group_ind == "I_0"]
  taxa_I1 <- tax_df$Pseudo_taxon[tax_df$group_ind == "I_1"]
  train_taxa <- c(taxa_I0, taxa_I1)
  train_pruned <- prune_taxa(train_taxa, train_final)
  
  return(list(
    pseudo_lfc = log_fold_change,
    Train = train_pruned
  ))
}

