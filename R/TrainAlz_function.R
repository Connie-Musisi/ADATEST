#' Train prediction model
#'
#' This function trains a linear regression model on microbiome training data to predict differential abundance.
#' It applies normalization to paired OTU counts, performs a sample-level permutation for the null group (`I_0`), 
#' and fits a regression model using sorted input features.
#' The resulting scores (`Train_parest`) are used downstream for test statistic computation.
#'
#' @param data A \code{phyloseq} object containing microbiome training data with taxonomic information and group labels.
#' @param group_var Character. The name of the grouping variable in the training metadata. Passed from \code{PseudoData()}.
#' @param group_levels Character vector of length two. The two levels of \code{group_var} to be compared. Passed from \code{PseudoData()}.
#' @param seed An optional seed-setting expression for reproducibility (default: \code{set.seed(19)}).
#'
#' @return A list containing:
#' \describe{
#'   \item{y}{A numeric vector of outcome labels: \code{-1} for \code{I_0} taxa, \code{1} for \code{I_1} taxa.}
#'   \item{Train_scaled}{A matrix of normalized and sorted training features.}
#'   \item{Train_results}{A data frame containing model predictions and the true outcome labels.}
#'   \item{Train_parest}{A numeric vector of trained regression coefficients (scores).}
#' }
#'
#' @details
#' This function is intended for internal use within the ADATEST pipeline. It constructs a feature matrix 
#' from paired taxa OTU counts using a normalization and sorting strategy. A simple linear model is fit 
#' to predict binary group membership, and the learned coefficients are interpreted as discriminative scores.
#'
#' @import phyloseq
#' @importFrom stats lm
#' @keywords internal
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you already constructed a pseudo-training object with group labels:
#' result <- Train_analyze(training_phyloseq_object)
#' print(result$Train_parest)
#' }

Train_analyze <- function(data,group_var,group_levels,seed){
  ########## Train data ######
  ## ---------------------------------------------------
  ##### Permutation #####
  # For I0 group
  train_I0 <- subset_taxa(data, group_ind == "I_0")
  otu_train_I0 <- otu_table(train_I0, taxa_are_rows = FALSE)
  tax_train_I0 <- tax_table(train_I0)

  # For I1 group
  train_I1 <- subset_taxa(data, group_ind == "I_1")
  
  # Convert OTU Table for I0 to DataFrame
  otu_train_I0 <- as.data.frame(otu_train_I0)
  
  # Permute the taxa labels in the OTU table for the I0 group for reproducibility
  Train_perm_I0<-otu_train_I0
  colnames(Train_perm_I0) <- taxa_names(train_I0) # not necessary
  rownames(Train_perm_I0) <- rownames(otu_train_I0)
  
  # Convert back to Phyloseq OTU Table format
  Train_otu_tab_I0 <- otu_table(Train_perm_I0, taxa_are_rows = FALSE)
  
  train_I0_perm<-phyloseq(Train_otu_tab_I0,sample_data(data),tax_train_I0)
  train_phy <- merge_phyloseq(train_I0_perm,train_I1)
  
  samdata <- as.data.frame(sample_data(train_phy))
  grp_vec <- samdata[[group_var]]
  n1 <- length(rownames(samdata)[grp_vec == group_levels[1]])
  n2 <- length(rownames(samdata)[grp_vec == group_levels[2]])

  
  
  ## ---------------------------------------------------------------------------------------------------------------------
  ############### Linear REGRESSION MODEL #################
  # 1. Extract data
  otu_abundance1 <- otu_table(train_phy) # OTU abundance matrix
  
  taxon_data_df1 <- data.frame(tax_table(train_phy))
  outcome_variable1 <- taxon_data_df1$group_ind 
  
  # 2. Prepare data for regression
  # Combine sample data with OTU abundance data
  regression_data1 <- as.matrix(t(otu_abundance1))
  test_func = function(x){regression_normalized<-Normalize(x[1:n1],
                                                           x[(n1+1):(n1+n2)])
  reg11 <- regression_normalized$sample1
  reg21 <- regression_normalized$sample2
  avg_reg11 <- mean(reg11)
  avg_reg21 <- mean(reg21)
  if(avg_reg21 < avg_reg11){
    return(c(reg21, reg11))
  }
  else {
    return(c(reg11, reg21))
  }}
  
  regression_data1 <- as.data.frame(t(apply(regression_data1,1,test_func)))
  
  regression_data1 = t(apply(regression_data1,1,function(x){c(sort(as.numeric(x[1:n1]),
                                                                   decreasing = FALSE),sort(as.numeric(x[c((n1+1):(n1+n2))]),
                                                                                            decreasing = FALSE))}))
  
  y <- numeric(nrow(regression_data1))
  y[outcome_variable1=="I_0"] = -1 # used to be zero
  y[outcome_variable1=="I_1"] = 1
  x <- as.matrix(regression_data1)
  
  if (any(is.na(y))) {
    stop("There are NA values in y")
  }
  if (length(y) != nrow(x)) {
    stop("The length of y2 does not match the number of rows in x1")
  }
  
  # Fitting linear model to train data
  Train_parest<-t(x)%*%y/length(y)
  Train_pred<-x%*%Train_parest
  
  ## ---------------------------------------------------------------------------------------------------------------------
  # Predicting the outcome
  
  Train_parest <- as.numeric(Train_parest)
  
  # matrix with both original group and predicted group
  colnames(Train_pred)[1] <- 'Predictions'
  Train_results <- cbind(Train_pred, y)
  
  # splitting the data into the 2 groups I_) and I_1 to get the distributions
  Train_results <- as.data.frame(Train_results)
  
  return(list(y=y,
              Train_scaled=regression_data1,
              Train_results=Train_results,
              Train_parest=Train_parest))
}
