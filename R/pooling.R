#' Fit Poolwise Regression Models
#'
#' Functions for calculating power and fitting regression models in studies
#' where a biomarker is measured in "pooled" samples rather than for each
#' individual.
#'
#' \tabular{ll}{
#' Package: \tab pooling \cr
#' Type: \tab Package \cr
#' Version: \tab 1.1.1 \cr
#' Date: \tab 2018-03-30 \cr
#' License: \tab GPL-3 \cr
#' }
#'
#' @author Dane R. Van Domelen \cr \email{vandomed@@gmail.com}
#'
#'
#' @references
#' Acknowledgment: This material is based upon work supported by the National
#' Science Foundation Graduate Research Fellowship under Grant No. DGE-0940903.
#'
#'
#' @docType package
#' @importFrom pracma hessian
#' @importFrom stats nlminb
#' @importFrom cubature adaptIntegrate
#' @importFrom mvtnorm dmvnorm
#' @import dvmisc
#' @name pooling
NULL
