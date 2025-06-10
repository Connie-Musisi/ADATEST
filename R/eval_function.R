#' Evaluate Performance Metrics for Differential Abundance Analysis
#'
#' This function computes key performance metrics, including False Discovery Rate (FDR), Sensitivity, Specificity, Type I Error Rate, and Power,
#' based on the classification of taxa as differentially abundant (DA) or not.
#'
#' @param Res A list containing results from the original dataset, including true and newly identified DA statuses.
#' @param isDA Logical vector. The ground truth for differential abundance status in the data (default: `Res$isDA`).
#' @param new_isDA Character vector. The newly classified differential abundance status (default: `Res$new_ind`).
#'
#' @return A list containing:
#' 	\item{FDR}{False Discovery Rate, computed as FP / (FP + TP).}
#' 	\item{Sensitivity}{True Positive Rate (TPR), computed as TP / (TP + FN).}
#' 	\item{Specificity}{True Negative Rate (TNR), computed as TN / (TN + FP).}
#' 	\item{Type_1_error}{Type I Error Rate, computed as FP / (FP + TN).}
#' 	\item{Power}{Statistical power, computed as 1 - beta, where beta is FN / (FN + TP).}
#'
#' @export
eval <- function(Res, isDA=Res$isDA, new_isDA=Res$new_ind){
  FP <- sum(isDA == "FALSE" & new_isDA == "I_1")
  TP <- sum(isDA == "TRUE" & new_isDA == "I_1")
  FN <- sum(isDA == "TRUE" & new_isDA == "I_0")
  TN <- sum(isDA == "FALSE" & new_isDA == "I_0")
  
  ## ----FDR------------------------------------------
  FDR <- FP/(FP+TP)
  
  ## ----Sensitivity----------------------------------
  TPR <- TP/(TP+FN)
  
  ## ----Specificity----------------------------------
  TNR <- TN/(TN+FP)
  
  ## ----Type 1 error rate----------------------------
  Type_1_error <- FP/(FP+TN)
  
  ## ----Power----------------------------------------
  beta <- FN/(FN+TP)
  
  return(list(FDR=FDR,
              Sensitivity=TPR,
              Specificity=TNR,
              Type_1_error=Type_1_error,
              Power=1-beta))
}
