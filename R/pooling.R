#' Fit Poolwise Regression Models
#'
#' Functions for calculating power and fitting regression models in studies
#' where a biomarker is measured in "pooled" samples rather than for each
#' individual.
#'
#' Currently, there are two functions for estimating the covariate-adjusted
#' log-odds ratio relating a continuous exposure measured in pools and subject
#' to errors and a binary outcome. The function \code{\link{p.logreg.xerrors}}
#' implements homogeneous-pools logistic regression, while
#' \code{\link{p.dfa.xerrors}} implements the discriminant function approach.
#'
"_PACKAGE"
