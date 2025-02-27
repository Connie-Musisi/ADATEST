#' Normalize two samples using percentile-based scaling/normalization
#'
#' This function normalizes two input samples by performing a percentile-based scaling. It pools the two samples together, 
#' calculates the 0th and 99th percentiles of the pooled data, and then scales the values within each sample. 
#' An additional adjustment is made based on the median values of the two samples.
#'
#' @param sample1 A numeric vector representing the first sample of data.
#' @param sample2 A numeric vector representing the second sample of data.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{sample1}{The normalized version of the first sample.}
#'   \item{sample2}{The normalized version of the second sample.}
#' }
#' @details
#' The function pools the two samples to calculate quantiles (0th and 99th) and uses the range between these 
#' quantiles as the normalization factor. If there are many zeros, the function defaults to using half the maximum value of the pooled data 
#' as the normalization factor. The normalized samples are then adjusted by the average of the medians of the two input samples.
#'
#' @examples
#' sample1 <- c(1, 2, 3, 4, 5)
#' sample2 <- c(2, 3, 4, 5, 6)
#' result <- Normalize(sample1, sample2)
#' result$sample1
#' result$sample2
#'
#' @export
Normalize<-function(sample1,sample2) {
  pooled<-c(sample1,sample2)
  q<-quantile(pooled,probs=c(0,0.99))
  norm.factor<-max((q[2]-q[1]),max(pooled)/2) # in case of many zeroes
  sample1_norm<-(sample1-q[1])/norm.factor
  sample2_norm<-(sample2-q[1])/norm.factor
  # version 2 needs following three lines
  #m<-(median(as.numeric(sample2))+median(as.numeric(sample1)))/2
  #sample1_norm<-sample1_norm-m
  #sample2_norm<-sample2_norm-m
  return(list(
    sample1=sample1_norm,
    sample2=sample2_norm
  ))
}
