#' Remove Taxa with Missing Values from a Phyloseq Object
#'
#' This function removes taxa (OTUs) that contain missing values (NA) from a given phyloseq object.
#' It ensures that only taxa with complete data are retained in both the OTU and taxonomy tables.
#'
#' @param physeq A phyloseq object containing OTU and taxonomy tables.
#'
#' @return A modified phyloseq object with taxa containing NAs removed.
#'
#' @export
NA.remove <- function(physeq) {
  otu <- otu_table(physeq)
  tax <- tax_table(physeq)
  
  # Check if taxa are rows or columns and remove NA accordingly
  if(taxa_are_rows(physeq)) {
    # Identify taxa (rows) that don't contain any NA
    valid_taxa <- !rowSums(is.na(otu))
    # Subset the OTU and taxa tables to keep only valid taxa
    otu <- otu[valid_taxa, ]
    tax <- tax[valid_taxa, ]
  } else {
    # Identify taxa (columns) without any NA
    valid_taxa <- !colSums(is.na(otu))
    # Subset the OTU and taxa tables to keep only valid taxa
    otu <- otu[, valid_taxa]
    tax <- tax[valid_taxa, ]
  }
  
  # Put the new tables back into the phyloseq object
  otu_table(physeq) <- otu
  tax_table(physeq) <- tax
  
  return(physeq)
}
