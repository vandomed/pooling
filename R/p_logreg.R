#' Poolwise Logistic Regression
#'
#' Fit homogeneous-pools logistic regression model described by Weinberg &
#' Umbach (1999).
#'
#'
#' @param g Numeric vector with pool sizes, i.e. number of members in each pool.
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
#' @param prev Numeric value specifying disease prevalence, allowing
#' for valid estimation of the intercept with case-control sampling. Can specify
#' \code{samp_y1y0} instead if sampling rates are known.
#'
#' @param samp_y1y0 Numeric vector of length 2 specifying sampling probabilities
#' for cases and controls, allowing for valid estimation of the intercept with
#' case-control sampling. Can specify \code{prev} instead if it's easier.
#'
#' @param estimate_var Logical value for whether to return variance-covariance
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
#' Weinberg, C.R. and Umbach, D.M. (1999) "Using pooled exposure assessment to
#' improve efficiency in case-control studies." \emph{Biometrics} \strong{55}:
#' 718--726.
#'
#' Weinberg, C.R. and Umbach, D.M. (2014) "Correction to 'Using pooled exposure
#' assessment to improve efficiency in case-control studies' by Clarice R.
#' Weinberg and David M. Umbach; 55, 718--726, September 1999."
#' \emph{Biometrics} \strong{70}: 1061.
#'
#'
#' @export
p_logreg <- function(g, y, x,
                     method = "glm", prev = NULL, samp_y1y0 = NULL,
                     estimate.var = TRUE, ...) {

  # Check that inputs are valid
  if (! method %in% c("glm", "ml")) {
    stop("The input 'method' should be set to 'glm' for glm function or 'ml' for maximum likelihood.")
  }
  if (! is.null(prev)) {
    if (prev < 0 | prev > 1) {
      stop("The input 'prev' is the disease prevalence, and must be between 0 and 1.")
    }
  }
  if (! is.null(samp_y1y0)) {
    if (! (length(samp_y1y0) == 2 & sum(samp_y1y0) == 1 &
           min(samp_y1y0) > 0 & max(samp_y1y0) < 1)) {
      stop("The input 'samp_y1y0' is the sampling probabilities for cases and controls, and should be a numeric vector of two probabilities adding to 1.")
    }
  }
  if (! is.logical(estimate.var)) {
    stop("The input 'estimate.var' should be TRUE or FALSE.")
  }

  # Get name of x input
  x.varname <- deparse(substitute(x))
  if (grep("$", x.varname)) {
    x.varname <- substr(x.varname,
                        start = which(unlist(strsplit(x.varname, "")) == "$") + 1,
                        stop = nchar(x.varname))
  }

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

  # Calculate offsets according to Weinberg and Umbach formula, incorporating
  # disease prevalence or sampling probabilities if known
  n <- length(y)
  locs.cases <- which(y == 1)
  n_1 <- sum(g[locs.cases])
  n_0 <- sum(g[-locs.cases])
  g.vals <- unique(g)
  qg <- rep(NA, n)

  if (! is.null(prev)) {

    for (jj in 1: length(g.vals)) {
      g.jj <- g.vals[jj]
      locs.g <- which(g == g.jj)
      n.casepools <- sum(g == g.jj & y == 1)
      n.controlpools <- sum(g == g.jj & y == 0)
      qg[locs.g] <- log(n.casepools / n.controlpools) -
        g.jj * log(prev / (1 - prev))
    }

  } else if (! is.null(samp_y1y0)) {

    for (jj in 1: length(g.vals)) {
      g.jj <- g.vals[jj]
      locs.g <- which(g == g.jj)
      n.casepools <- sum(g == g.jj & y == 1)
      n.controlpools <- sum(g == g.jj & y == 0)
      qg[locs.g] <- log(n.casepools / n.controlpools) -
        g.jj * log(n_1 / n_0) - g.jj * log(samp_y1y0[2] / samp_y1y0[1])
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
      extra.args$start <- rep(0.01, n.betas)
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
