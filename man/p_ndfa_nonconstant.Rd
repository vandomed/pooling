% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p_ndfa_nonconstant.R
\name{p_ndfa_nonconstant}
\alias{p_ndfa_nonconstant}
\title{Normal Discriminant Function Approach for Estimating Odds Ratio with Exposure
Measured in Pools and Potentially Subject to Additive Normal Errors
(Non-constant Odds Ratio Version)}
\usage{
p_ndfa_nonconstant(g, y, xtilde, c = NULL, errors = "processing",
  start_nonvar_var = c(0.01, 1), lower_nonvar_var = c(-Inf, 1e-04),
  upper_nonvar_var = c(Inf, Inf), jitter_start = 0.01,
  nlminb_list = list(control = list(trace = 1, eval.max = 500, iter.max =
  500)), hessian_list = list(method.args = list(r = 4)),
  nlminb_object = NULL)
}
\arguments{
\item{g}{Numeric vector of pool sizes, i.e. number of members in each pool.}

\item{y}{Numeric vector of poolwise Y values (number of cases in each pool).}

\item{xtilde}{Numeric vector (or list of numeric vectors, if some pools have
replicates) with Xtilde values.}

\item{c}{Numeric matrix with poolwise \strong{C} values (if any), with one
row for each pool. Can be a vector if there is only 1 covariate.}

\item{errors}{Character string specifying the errors that X is subject to.
Choices are \code{"neither"}, \code{"processing"} for processing error only,
\code{"measurement"} for measurement error only, and \code{"both"}.}

\item{start_nonvar_var}{Numeric vector of length 2 specifying starting value
for non-variance terms and variance terms, respectively.}

\item{lower_nonvar_var}{Numeric vector of length 2 specifying lower bound for
non-variance terms and variance terms, respectively.}

\item{upper_nonvar_var}{Numeric vector of length 2 specifying upper bound for
non-variance terms and variance terms, respectively.}

\item{jitter_start}{Numeric value specifying standard deviation for mean-0
normal jitters to add to starting values for a second try at maximizing the
log-likelihood, should the initial call to \code{\link[stats]{nlminb}} result
in non-convergence. Set to \code{NULL} for no second try.}

\item{nlminb_list}{List of arguments to pass to \code{\link[stats]{nlminb}}
for log-likelihood maximization.}

\item{hessian_list}{List of arguments to pass to
\code{\link[numDeriv]{hessian}} for approximating the Hessian matrix. Only
used if \code{estimate_var = TRUE}.}

\item{nlminb_object}{Object returned from \code{\link[stats]{nlminb}} in a
prior call. Useful for bypassing log-likelihood maximization if you just want
to re-estimate the Hessian matrix with different options.}
}
\value{
List containing:
\enumerate{
\item Numeric vector of parameter estimates.
\item Variance-covariance matrix.
\item Returned \code{\link[stats]{nlminb}} object from maximizing the
log-likelihood function.
\item Akaike information criterion (AIC).
}
}
\description{
Assumes exposure given covariates and outcome is a normal-errors linear
regression. Pooled exposure measurements can be assumed precise or subject to
additive normal processing error and/or measurement error. Parameters are
estimated using maximum likelihood.
}
\references{
Lyles, R.H., Van Domelen, D.R., Mitchell, E.M. and Schisterman, E.F. (2015)
"A discriminant function approach to adjust for processing and measurement
error When a biomarker is assayed in pooled samples."
\emph{Int. J. Environ. Res. Public Health} \strong{12}(11): 14723--14740.

Schisterman, E.F., Vexler, A., Mumford, S.L. and Perkins, N.J. (2010) "Hybrid
pooled-unpooled design for cost-efficient measurement of biomarkers."
\emph{Stat. Med.} \strong{29}(5): 597--613.
}
