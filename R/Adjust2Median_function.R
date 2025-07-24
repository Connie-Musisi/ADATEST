#' Adjust Counts to Median Fold-Change Between Two Groups
#'
#' This function scales the counts of samples in the reference group (group == 0)
#' so that the overall abundance shift between the two groups matches the median
#' log2-fold change observed across taxa. This adjustment can help mitigate
#' global compositional effects when comparing differential abundance.
#'
#' @param data A \code{phyloseq} object. It must contain an OTU table
#'   accessible via \code{otu_table(data)} and sample metadata including a
#'   \code{group} factor or numeric vector with values 0 and 1 (two groups).
#'
#' @return A \code{phyloseq} object with the same taxa and sample metadata,
#'   but with the OTU counts for samples in \code{group == 0} multiplied by
#'   the median fold-change factor.
#'
#' @details
#' The function works as follows:
#' 1. Transpose the OTU count matrix so that taxa become columns and samples rows.
#' 2. Compute, for each taxon, the smallest non-zero count across all samples.
#' 3. Calculate log2 fold-changes (LFC) between the two groups (1 vs. 0)
#'    by adding the pseudocounts to avoid zeros:
#'    \code{lfc = log2((mean(group==1) + min.nonzero) / (mean(group==0) + min.nonzero))}.
#' 4. Compute the median of the back-transformed fold-changes: \code{med0 = 2^(median(lfc))}.
#' 5. Multiply all counts in group 0 by \code{med0} to align the median abundance shift.
#' 6. Reconstruct a new \code{phyloseq} object with adjusted OTU counts.
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' # Assume `physeq` is a phyloseq object with sample_data containing `group` (0/1)
#' physeq_adj <- Adjust2Median(physeq)
#' # Inspect adjusted counts
#' otu_table(physeq_adj)
#' }
#'
#' @export
Adjust2Median <- function(data,group_var,group_levels) {
  # Extract and transpose OTU count matrix
  if(taxa_are_rows(data)){
    tmp <- t(otu_table(data))
  }
  # Extract grouping variable (0/1)
  samdata <- as.data.frame(sample_data(data))
  grp_vec <- samdata[[group_var]]
  
  # Compute pseudocount for each taxon as the minimum non-zero count
  min.without.zero <- apply(tmp, 2, FUN = function(x) {
    min(x[x != 0], na.rm = TRUE)
  })
  
  # Calculate log2 fold-changes between groups, adding pseudocount
  lfc2 <- log2((colMeans(tmp[grp_vec == group_levels[2], ] + min.without.zero)) /
                 (colMeans(tmp[grp_vec == group_levels[1], ] + min.without.zero)))
  
  # Calculate median fold-change factor
  med0 <- 2^median(lfc2)
  
  # Scale counts in group 0 by median fold-change
  tmp[grp_vec == group_levels[1], ] <- tmp[grp_vec == group_levels[1], ] * med0
  
  # Reconstruct phyloseq object with adjusted OTU counts
  physeq <- phyloseq(t(tmp), tax_table(data), sample_data(data))
  
  return(physeq)
}
