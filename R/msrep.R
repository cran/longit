#' @title longitudinal data
#'
#' @description Longitudinal observation on single variable at different timepoints. Observations arranged in a column as the patient with corresponding column of ID.
#' @usage data(msrep)
#' @format A \code{tibble} with 7 columns which are :
#' \describe{
#' \item{Subject}{Patient ID}
#' \item{Gender}{Categorical numeric variable, 1 if Males and 0 if female}
#' \item{Age}{Time or age at which observations were taken from every subjects}
#' \item{x1,...,x4}{Columns stating number of observations at age 18,10,12 and 14}}
#'
"msrep"
