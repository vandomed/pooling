#' Poolwise Logistic Regression
#'
#' Fit logistic regression model with variables measured in pools rather than
#' individual participants. Requires that each pool is comprised of all cases or
#' all controls.
#'
#' \emph{Statistical setup:}
#'
#' The goal is to estimate coefficients in an individual-level logistic
#' regression model relating a continuous outcome \code{Y} to covariates
#' \code{\strong{X}}:
#'
#' \code{logit[P(Y_ij = 1)] = beta_0 + \bold{beta_x}^T \bold{X_ij}}
#'
#' But we observe poolwise values rather than individual values. In the
#' \code{i}th pool, which has \code{g_i} members that are either all cases or
#' all controls, we observe the following variables:
#'
#' \code{Y_i = 1} if all \code{Y_ij = 1} and \code{0} if all \code{Y_ij = 0,
#' j = 1, ..., g_i} \cr
#' \code{\bold{X_i} = sum_{j=1}^{g_i} \bold{X_ij}}
#'
#' A corresponding poolwise model has been described by Weinberg \& Umbach
#' (1999, 2014) and Lyles et al. (2016):
#'
#' \code{logit[P(Y_i = 1)] = g_i beta_0 + \bold{beta_x}^T \bold{X_i} + qg_i}
#'
#' Here \code{qg_i} is an offset defined as:
#'
#' \code{qg_i = ln(\# case pools of size g_i) / \# control pools of size g_i) +
#' g_i ln(n_0 / n_1) + ln[P(A|Y=1) / P(A|Y=0)]}
#'
#' where \code{n_0 (n_1)} represent the number of individual controls (cases) in
#' the study, across all pools, and \code{P(A|Y=1)} and \code{P(A|Y=0)} are
#' accrual probabilities for cases and controls. With prospective sampling,
#' these latter probabilities are equal and so the third term in the offset
#' drops out. With case-control sampling, these probabilities can be specified
#' through the \code{p_sample.y1y0} input, or alternatively (and equivalently)
#' the disease prevalence can be specified through the \code{p_y1} input. If
#' neither input is specified, the \code{beta_0} estimate will not be valid.
#'
#' Thus, with poolwise data, one can fit a no-intercept logistic regression
#' model with \code{g_i} included as a covariate and with an offset of
#' \code{qg_i}. If \code{method = "glm"}, this approach is taken.
#'
#' If \code{method = "ml"}, maximum likelihood is used instead. The likelihood
#' contribution for each observed poolwise \code{(Y_i, \bold{X_i})} is
#' \code{f(Y_i, \bold{X_i})}, which is proportional to \code{f(Y_i|\bold{X_i})},
#' which is the Bernoulli probability mass function corresponding to
#' \code{Y_i|\bold{X_i} ~ Bern(p_i = [g_i beta_0* + \bold{beta_x}^T \bold{X_i} +
#' ln(qg_i)]^(-1))}.
#'
#'
#' @param g Numeric vector with pool sizes (number of members in each pool).
#'
#' @param y Numeric vector with poolwise \code{Y} values, coded 0 if all members
#' are controls and 1 if all members are cases.
#'
#' @param x Numeric matrix with poolwise \strong{\code{X}} values, with one row
#' for each pool. Can be a vector if there is only 1 predictor.
#'
#' @param method Character string specifying method to use for estimation.
#' Choices are "glm" for \code{\link[stats]{glm}} function and \code{"ml"} for
#' maximum likelihood.
#'
#' @param p_y1 Numeric value specifying disease prevalence, allowing
#' for valid estimation of the intercept with case-control sampling. Can specify
#' \code{p_sample.y1y0} instead if sampling rates are known.
#'
#' @param p_sample.y1y0 Numeric vector of length 2 specifying sampling
#' probabilities for cases and controls, allowing for valid estimation of the
#' intercept with case-control sampling. Can specify \code{p_y1} instead if it's
#' easier.
#'
#' @param estimate.var Logical value for whether to return variance-covariance
#' matrix for parameter estimates.
#'
#' @param ... Additional arguments to pass to \code{\link[stats]{nlminb}}.
#'
#'
#' @return
#' List containing:
#' \enumerate{
#' \item Numeric vector of parameter estimates.
#' \item Variance-covariance matrix (if \code{estimate.var = TRUE}).
#' \item Fitted \code{\link[stats]{glm}} object (if \code{method = "glm"}) or
#' returned \code{\link[stats]{nlminb}} object (if \code{method = "ml"}).
#' \item Akaike information criterion (AIC).
#' }
#'
#'
#' @references
#' Lyles, R.H., Mitchell, E.M., Weinberg, C.R., Umbach, D.M. and Schisterman,
#' E.F. (2016) "An Efficient Design Strategy for Logistic Regression Using
#' Outcome-Dependent Pooling of Biospecimens Prior to Assay." \emph{Biometrics}
#' \strong{72}: 965--975.
#'
#' Weinberg, C.R. and Umbach, D.M. (1999) "Using Pooled Exposure Assessment to
#' Improve Efficiency in Case-Control Studies." \emph{Biometrics} \strong{55}:
#' 718--726.
#'
#' Weinberg, C.R. and Umbach, D.M. (2014) "Correction to 'Using Pooled Exposure
#' Assessment to Improve Efficiency in Case-Control Studies' by Clarice R.
#' Weinberg and David M. Umbach; 55, 718--726, September 1999."
#' \emph{Biometrics} \strong{70}: 1061.
#'
#' Acknowledgment: This material is based upon work supported by the National
#' Science Foundation Graduate Research Fellowship under Grant No. DGE-0940903.
#'
#'
#' @export
p_logreg <- function(g, y, x,
                     method = "glm", p_y1 = NULL, p_sample.y1y0 = NULL,
                     estimate.var = TRUE, ...) {

  # Check that inputs are valid
  if (! method %in% c("glm", "ml")) {
    stop("The input 'method' should be set to 'glm' for glm function or 'ml' for maximum likelihood.")
  }
  if (! is.null(p_y1)) {
    if (p_y1 < 0 | p_y1 > 1) {
      stop("The input 'p_y1' is the disease prevalence, and must be between 0 and 1.")
    }
  }
  if (! is.null(p_sample.y1y0)) {
    if (! (length(p_sample.y1y0) == 2 & sum(p_sample.y1y0) == 1 &
           min(p.sample.y1y0) > 0 & max(p.sample.y1y0) < 1)) {
      stop("The input 'p_sample_y1y0' is the sampling probabilities for cases and controls, and should be a numeric vector of two probabilities adding to 1.")
    }
  }
  if (! is.logical(estimate.var)) {
    stop("The input 'estimate.var' should be TRUE or FALSE.")
  }

  # Get name of x input
  x.varname <- deparse(substitute(x))

  # Get number of X variables (and assign names)
  if (class(x) != "matrix") {
    x <- as.matrix(x)
  }
  n.xvars <- ncol(x)
  x.varnames <- colnames(x)
  if (is.null(x.varnames)) {
    if (n.xvars == 1) {
      x.varnames <- x.varname
    } else {
      x.varnames <- paste("x", 1: n.xvars, sep = "")
    }
  }

  # Create matrix of (g, X) values
  gx <- cbind(g, x)

  # Calculate offsets according to Weinberg and Umbach (Biometrics 1999, 2010)
  # formula, incorporating disease prevalence or sampling probabilities if known
  n.pools <- length(y)
  locs.cases <- which(y == 1)
  n_1 <- sum(g[locs.cases])
  n_0 <- sum(g[-locs.cases])
  g.vals <- unique(g)
  qg <- rep(NA, n.pools)

  if (! is.null(p_y1)) {

    for (jj in 1: length(g.vals)) {
      g.jj <- g.vals[jj]
      locs.g <- which(g == g.jj)
      n.casepools <- sum(g == g.jj & y == 1)
      n.controlpools <- sum(g == g.jj & y == 0)
      qg[locs.g] <- log(n.casepools / n.controlpools) -
        g.jj * log(p_y1 / (1 - p_y1))
    }

  } else if (! is.null(p_sample.y1y0)) {

    for (jj in 1: length(g.vals)) {
      g.jj <- g.vals[jj]
      locs.g <- which(g == g.jj)
      n.casepools <- sum(g == g.jj & y == 1)
      n.controlpools <- sum(g == g.jj & y == 0)
      qg[locs.g] <- log(n.casepools / n.controlpools) -
        g.jj * log(n_1 / n_0) - g.jj * log(p_sample.y1y0[2] /
                                             p_sample.y1y0[1])
    }

  } else {

    for (jj in 1: length(g.vals)) {
      g.jj <- g.vals[jj]
      locs.g <- which(g == g.jj)
      n.casepools <- sum(g == g.jj & y == 1)
      n.controlpools <- sum(g == g.jj & y == 0)
      qg[locs.g] <- log(n.casepools / n.controlpools) - g.jj * log(n_1 / n_0)
    }

  }

  # Create labels for parameter estimates
  beta.labels <- paste("beta", c("0", x.varnames), sep = "_")

  # Use glm function or ML depending on method input
  if (method == "glm") {

    # Fit logistic regression
    glm.fit <- glm(y ~ gx - 1, family = "binomial", offset = qg)

    # Create list to return
    theta.hat <- glm.fit$coef
    names(theta.hat) <- beta.labels
    ret.list <- list(theta.hat = theta.hat)

    # If requested, add variance-covariance matrix to ret.list
    if (estimate.var) {
      glm.variance <- vcov(glm.fit)
      colnames(glm.variance) <- rownames(glm.variance) <- beta.labels
      ret.list$glm.var <- glm.variance
    }

    # Add fitted glm object and AIC to ret.list
    ret.list$glm.fit <- glm.fit
    ret.list$aic <- AIC(glm.fit)

  } else if (method == "ml") {

    # Log-likelihood function
    n.betas <- length(beta.labels)
    ll.f <- function(f.theta) {

      # Likelihood:
      # L_i = f(Y|X)

      # P(Y|X)
      eta <- gx %*% f.theta + qg
      p_y.x <- (1 + exp(-eta))^(-1)

      # Log-likelihood
      ll <- sum(dbinom(x = y, size = 1, prob = p_y.x, log = TRUE))
      return(-ll)

    }

    # Create list of extra arguments, and assign default starting values and
    # lower values if not specified by user
    extra.args <- list(...)
    if (is.null(extra.args$start)) {
      extra.args$start <- rep(0, n.betas)
    }
    if (is.null(extra.args$lower)) {
      extra.args$lower <- rep(-Inf, n.betas)
    }
    if (is.null(extra.args$control$rel.tol)) {
      extra.args$control$rel.tol <- 1e-6
    }

    # Obtain ML estimates
    ml.max <- do.call(nlminb, c(list(objective = ll.f), extra.args))

    # Create list to return
    theta.hat <- ml.max$par
    names(theta.hat) <- beta.labels
    ret.list <- list(theta.hat = theta.hat)

    # If requested, add variance-covariance matrix to ret.list
    if (estimate.var) {
      hessian.mat <- pracma::hessian(f = ll.f, x0 = theta.hat)
      theta.variance <- solve(hessian.mat)
      colnames(theta.variance) <- rownames(theta.variance) <- beta.labels
      ret.list$theta.var <- theta.variance
    }

    # Add nlminb object and AIC to ret.list
    ret.list$nlminb.object <- ml.max
    ret.list$aic <- 2 * (length(theta.hat) + ml.max$objective)

  }

  # Return ret.list
  return(ret.list)

}
