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
