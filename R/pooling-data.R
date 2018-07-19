#' Dataset for Examples in p_linreg_yerrors
#'
#' List containing (1) data frame with poolwise (Y, X1, X2) values and (2) list
#' with replicate Y values. Data were generated under linear regression of Y on
#' (X1, X2) with beta_0 = 0.25, beta_x1 = 0.5, and beta_x2 = 0.25, with Y
#' subject to processing error and measurement error.
#'
#' @name dat_p_linreg_yerrors
#' @docType data
#' @source Simulated data in R.
NULL
#'
#'
#' Dataset with Simulated (Y, C) Values for Examples in dfa_xerrors and
#' logreg_xerrors
#'
#' Includes 5,000 observations. The \code{Xtilde} values are stored as a
#' separate list called \code{\link{dat1_xtilde}}. The data were generated with
#' a true log-odds ratio of 0.5 for \code{X} and \code{Y}, adjusted for
#' \code{C}. The \code{Xtilde} measurements are subject to measurement error.
#'
#' @name dat1
#' @docType data
#' @source Simulated data in R.
NULL
#'
#'
#' Dataset with Simulated Xtilde Values for Examples in dfa_xerrors and
#' logreg_xerrors
#'
#' Includes 5,000 observations, 30 of which have replicates. The \code{(Y, C)}
#' values are stored as a separate data frame called \code{\link{dat1}}. The
#' data were generated with a true log-odds ratio of 0.5 for \code{X} and
#' \code{Y}, adjusted for \code{C}. The \code{Xtilde} measurements are subject
#' to measurement error.
#'
#' @name dat1_xtilde
#' @docType data
#' @source Simulated data in R.
NULL
#'
#'
#' Dataset with Simulated (Y, Xtilde, C) Values for Examples in p_dfa_xerrors
#' and p_logreg_xerrors
#'
#' Includes 4,999 pooled observations, with a roughly equal number of pools of
#' size 1, 2, and 3. The data were generated with a true log-odds ratio of 0.5
#' for \code{X} and \code{Y}, adjusted for \code{C}. The \code{Xtilde}
#' measurements are subject to processing error.
#'
#' @name pdat1
#' @docType data
#' @source Simulated data in R.
NULL
#'
#'
#' Dataset with Simulated (Y, Xtilde) Values for Examples in p_dfa_xerrors2 and
#' p_logreg_xerrors2
#'
#' Includes 248 pooled observations, with a roughly equal number of pools of
#' size 1, 2, and 3. The individual-level \code{C} values are stored as a
#' separate list called \code{\link{pdat2_c}}. The data were generated with
#' a true log-odds ratio of 0.5 for \code{X} and \code{Y}, adjusted for
#' \code{C}. The \code{Xtilde} measurements are subject to processing error.
#'
#' @name pdat2
#' @docType data
#' @source Simulated data in R.
NULL
#'
#'
#' Dataset with Simulated C Values for Examples in p_dfa_xerrors2 and
#' p_logreg_xerrors2
#'
#' Includes 248 sets of individual-level \code{C} values. The \code{(Y, Xtilde)}
#' values are stored as a separate data frame called \code{\link{pdat2}}. The
#' data were generated with a true log-odds ratio of 0.5 for \code{X} and
#' \code{Y}, adjusted for \code{C}. The \code{Xtilde} measurements are subject
#' to processing error.
#'
#' @name pdat2_c
#' @docType data
#' @source Simulated data in R.
NULL
#'
#'
#' Dataset with Simulated Values for Examples in cond_logreg
#'
#' Includes 50 pools of size 1 and 75 pools of size 2. The error-prone
#' \code{Xtilde} values for case pools and control pools are stored as separate
#' lists called \code{\link{xtilde1_matched}} and \code{\link{xtilde0_matched}}.
#'
#' @name pdat_matched
#' @docType data
#' @source Simulated data in R.
NULL
#'
#'
#' Dataset with Simulated Xtilde1's for Examples in cond_logreg
#'
#' Includes 50 case pools of size 1 (25 of which have replicates) and 75 of
#' size 2. Other variables are stored in a separate data frame called
#' \code{\link{pdat_matched}}.
#'
#' @name xtilde1_matched
#' @docType data
#' @source Simulated data in R.
NULL
#'
#'
#' Dataset with Simulated Xtilde0's for Examples in cond_logreg
#'
#' Includes 50 control pools of size 1 (25 of which have replicates) and 75 of
#' size 2. Other variables are stored in a separate data frame called
#' \code{\link{pdat_matched}}.
#'
#' @name xtilde0_matched
#' @docType data
#' @source Simulated data in R.
NULL
