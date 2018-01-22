#' @title
#' Checks whether true value is between CI limits.
#'
#' @description
#' This function takes a true value and checks whether the value is between
#' the upper and lower limits of a given Confidence Interval.
#'
#' @export
indicating <- function(UPPER, LOWER, trueVal){
  IND <- trueVal >= LOWER & trueVal <= UPPER
  COVERAGE <- apply(IND, 1, mean, na.rm=TRUE)
  return(COVERAGE)
}
