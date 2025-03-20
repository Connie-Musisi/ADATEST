#' Train prediction model
#'
#' This function trains a linear regression model on microbiome data to predict differential abundance in the training dataset.
#' It applies a permutation approach for taxa labels and normalizes OTU counts before model fitting.
#' The key output is `Train_parest`, a vector of trained regression parameters/partial least square estimates
#' used for prediction and we call these the scores.
#'
#' @param data A `phyloseq` object containing microbiome data with taxonomic and sample information.
#' @param seed An optional seed for reproducibility. Default is `set.seed(19)`.
#'
#' @return A list containing:
#' \describe{
#'   \item{y}{A numeric vector of outcome variables (`-1` for `I_0`, `1` for `I_1`).}
#'   \item{Train_scaled}{A matrix of scaled training data used for fitting the model.}
#'   \item{Train_results}{A data frame containing model predictions and actual outcomes.}
#'   \item{Train_parest}{A numeric vector of trained regression parameters.}
#' }
#'
#' @import phyloseq
#' @importFrom stats lm
#' @export
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data("GlobalPatterns") # Example dataset
#' result <- Train_analyze(GlobalPatterns)
#' print(result$Train_results)
#' }
Train_analyze <- function(data, seed=set.seed(19)){
  ########## Train data ######
  ## ---------------------------------------------------
  ##### Permutation #####
  # For I0 group
  train_I0 <- subset_taxa(data, group_ind == "I_0")
  otu_train_I0 <- otu_table(train_I0, taxa_are_rows = FALSE)
  tax_train_I0 <- tax_table(train_I0)
  n0<-nrow(otu_train_I0)/2 
  
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
  
  ## ---------------------------------------------------------------------------------------------------------------------
  ############### Linear REGRESSION MODEL #################
  # 1. Extract data
  otu_abundance1 <- otu_table(train_phy) # OTU abundance matrix
  
  taxon_data_df1 <- data.frame(tax_table(train_phy))
  outcome_variable1 <- taxon_data_df1$group_ind 
  
  # 2. Prepare data for regression
  # Combine sample data with OTU abundance data
  regression_data1 <- as.matrix(t(otu_abundance1))
  test_func = function(x){regression_normalized<-Normalize(x[1:n0],
                                                           x[(n0+1):(2*n0)])
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
  
  regression_data1 = t(apply(regression_data1,1,function(x){c(sort(as.numeric(x[1:n0]),
                                                                   decreasing = FALSE),sort(as.numeric(x[c((n0+1):(2*n0))]),
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
