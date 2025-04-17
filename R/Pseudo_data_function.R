#' PseudoData Function
#' 
#' This function constructs a pseudo-dataset from a given phyloseq object by generating 
#' artificial pseudo-taxa from the group with the largest sample size. It is designed to 
#' support differential abundance (DA) analysis while maintaining valid inference without bias.
#' 
#' @param simdata_filter A phyloseq object containing the microbiome data.
#' @param t0 Lower threshold for log fold change (LFC) classification of pseudo taxa (default: 1).
#' @param t1 Upper threshold for DA classification (default: 1.5).
#' @param t2 Maximum threshold for DA classification (default: +Inf).
#' @param n.taxa0 Number of taxa in the non-DA group, if known (if default (NULL), it is set dynamically 
#' based on half the total number of taxa.).
#' @param empirical_adjust Logical; if TRUE, adjusts the pseudo-dataset based on empirical distribution.
#' 
#' @return A list containing:
#' \item{pseudo_lfc}{Log fold change values of pseudo-taxa.}
#' \item{Train}{Training dataset for downstream analysis.}
#' 
#' @details
#' 1. **Pseudo-Dataset Construction**: The function selects the group with the largest 
#'    sample size and constructs pseudo-taxa by pairing taxa within this group.
#' 2. **Empirical Adjustment**: If `empirical_adjust = TRUE`, it downsamples data to match 
#'    the original group distribution.
#' 3. **Differential Abundance Classification**: Log fold changes (LFCs) are computed, and 
#'    pseudo-taxa are categorized into three subsets based on thresholds `t0`, `t1`, and `t2`:
#'    - `I_0`: Non-differentially abundant pseudo-taxa.
#'    - `I_1`: Differentially abundant pseudo-taxa.
#'    - `I_B`: Biologically irrelevant or extreme DA pseudo-taxa.
#' 4. **Training Dataset Generation**: The final dataset is created by normalizing and 
#'    partitioning pseudo-taxa based on LFC.
#' 
#' This function is useful for developing DA testing procedures that control false 
#' discovery rates while leveraging synthetic data for statistical power.
#' 
#' @examples
#' result <- PseudoData(simdata_filter)
#' pseudo_lfc <- result$pseudo_lfc
#' training_set <- result$Train
#' 
#' @export
PseudoData <- function(simdata_filter,
                       group_var,
                       group_levels,
                       t0=1, 
                       t1=1.5, 
                       t2=+Inf,
                       empirical_adjust = TRUE,
                       n.taxa0 = NULL) {
  
  
  samdata <- as.data.frame(sample_data(simdata_filter))
  grp_vec <- samdata[[group_var]]
  
  # Subset based on the group variable directly for group 1
  samps_grp1 <- rownames(samdata)[grp_vec == group_levels[1]]
  simdata_group1 <- prune_samples(samps_grp1, simdata_filter)
  if (taxa_are_rows(simdata_group1)) {
    otutab_group1 <- t(otu_table(simdata_group1))
  } else {
    otutab_group1 <- otu_table(simdata_group1)
  }
  n_samples_group1 <- nrow(otutab_group1)
  
  # Subset for group 2 using the same direct approach
  samps_grp2 <- rownames(samdata)[grp_vec == group_levels[2]]
  simdata_group2 <- prune_samples(samps_grp2, simdata_filter)
  if (taxa_are_rows(simdata_group2)) {
    otutab_group2 <- t(otu_table(simdata_group2))
  } else {
    otutab_group2 <- otu_table(simdata_group2)
  }
  n_samples_group2 <- nrow(otutab_group2)
  
  
  n.taxa <- ncol(otutab_group1)
  if (is.null(n.taxa0)) {
    n.taxa0 <- round(n.taxa / 2)
  }
  
  # Get the taxon names from the OTU table
  taxon_names <- taxa_names(simdata_group1)
  taxon_name_pairs <- c()
  
  use <- nrow(otutab_group1)
  for (i in 1:(n.taxa - 1)) {
    pairotu_subgroup1 <- t(matrix(c(otutab_group1[, i]), nrow = n.taxa - i, ncol = use, byrow = TRUE))
    pairotu_subgroup2 <- otutab_group1[, (i + 1):n.taxa]
    # Create the taxon name pairs
    taxon_name_pairs <- c(taxon_name_pairs, paste(rep(taxon_names[i], n.taxa - i), taxon_names[(i + 1):n.taxa], sep = "-"))
    
    if (i == 1) {
      combined_paired_otu <- rbind(pairotu_subgroup1, pairotu_subgroup2)
    } else {
      combined_paired_otu <- cbind(combined_paired_otu, rbind(pairotu_subgroup1, pairotu_subgroup2))
    }
  }
  
  # Remove taxa with too many zeroes 
  idx.keep.freq0 <- apply(combined_paired_otu, 2, function(x) { mean(x == 0) }) < 0.8
  idx.keep.avg <- apply(combined_paired_otu, 2, function(x) { mean(x) }) > 0.001
  idx.keep <- idx.keep.avg & idx.keep.freq0
  combined_paired_otu <- combined_paired_otu[, idx.keep]
  
  # Construct a New Phyloseq Object
  tax_data <- tax_table(as.matrix(taxon_name_pairs[idx.keep]))
  rownames(tax_data) <- taxon_name_pairs[idx.keep]
  colnames(tax_data) <- c("Pseudo_taxon")
  
  # Duplicate the sample data for the new structure
  samdata_new <- sample_data(simdata_group1)
  samdata_new <- as.data.frame(rbind(samdata_new, samdata_new))
  samdata_new <- sample_data(samdata_new)
  
  # Create a new OTU table with the combined paired OTU data (include names)
  otu_table_new <- otu_table(combined_paired_otu, taxa_are_rows = FALSE)
  colnames(otu_table_new) <- taxon_name_pairs[idx.keep]
  row.names(otu_table_new) <- row.names(samdata_new)
  
  training_dataset_ps <- phyloseq(otu_table_new, tax_data, samdata_new)
  training_data <- phyloseq(otu_table_new, tax_data, samdata_new)
  
  # Using empirical distribution to downsize group 2 training data to original
  if (empirical_adjust) {
    train_group2 <- otu_table_new[(n_samples_group1 + 1):(2 * n_samples_group1), ]
    
    redtrain_group2 <- matrix(nrow = n_samples_group2, ncol = ncol(train_group2))
    apply_emp <- function(x) {
      set.seed(19)
      emp_dist <- ecdf(x)
      quantile(emp_dist, probs = runif(n_samples_group2))
    }
    redtrain_group2 <- apply(train_group2, 2, apply_emp)
    
    train_otu_table <- otu_table(rbind(otu_table_new[1:n_samples_group1, ],
                                       redtrain_group2), taxa_are_rows = FALSE)
    rownames(train_otu_table) <- rownames(simdata_filter@sam_data)
    train_sample_data <- sample_data(simdata_filter)[, "group"]
    training_data <- phyloseq(train_otu_table, tax_data, train_sample_data)
  }
  
  # Log fold change calculation
  otu_train_filter <- training_data@otu_table
  data_group1 <- otu_train_filter[1:n_samples_group1, ]
  data_group2 <- otu_train_filter[(n_samples_group1 + 1):(n_samples_group1 + n_samples_group2), ]
  
  mean_group1 <- colMeans(data_group1)
  mean_group2 <- colMeans(data_group2)
  
  min.without.zero <- apply(otu_train_filter, 2, function(x) { min(x[x != 0]) })
  log_fold_change <- log2((mean_group2 + min.without.zero) / (mean_group1 + min.without.zero))
  # Calculate δ_med as the median of the log fold changes
  delta_med <- median(log_fold_change, na.rm = TRUE)
  
  # Center the log fold changes
  lfc_centered <- log_fold_change - delta_med
  
  # Extract the taxonomic table as a data frame
  tax_df <- as.data.frame(tax_table(training_data))
  tax_df$lfc <- log_fold_change
  tax_df$lfc_centered <- lfc_centered
  
  # Assign group indicators based on LFC using thresholds t0, t1, t2
  tax_df$group_ind <- ifelse(abs(lfc_centered) < t0, "I_0",
                             ifelse(abs(lfc_centered) > t1 & abs(lfc_centered) < t2, "I_1", "I_B"))
  
  sample_df <- as.data.frame(sample_data(training_data))
  sample_df$new_group <- c(rep(0, n_samples_group1), rep(1, n_samples_group2))
  sample_data_new <- sample_data(sample_df)
  
  train_final <- phyloseq(otu_table(training_data), sample_data_new, tax_table(as.matrix(tax_df)))
  
  taxa_I0 <- tax_df$Pseudo_taxon[tax_df$group_ind == "I_0"]
  taxa_I1 <- tax_df$Pseudo_taxon[tax_df$group_ind == "I_1"]
  train_taxa <- c(taxa_I0, taxa_I1)
  train_pruned <- prune_taxa(train_taxa, train_final)
  
  return(list(pseudo_lfc = log_fold_change,
              delta_med = delta_med,
              lfc_centered = lfc_centered,
              Train = train_pruned))
}
