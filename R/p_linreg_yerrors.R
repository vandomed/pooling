#' Linear Regression of Y vs. Covariates with Y Measured in Pools and
#' (Potentially) Subject to Additive Normal Errors
#'
#' Assumes outcome given covariates is a normal-errors linear regression. Pooled
#' outcome measurements can be assumed precise or subject to additive normal
#' processing error and/or measurement error. Replicates are supported.
#'
#' The individual-level model of interest for Y|\strong{X} is:
#'
#' Y = beta_0 + \strong{beta_x}^T \strong{X} + e, e ~ N(0, sigsq)
#'
#' The implied model for summed Y*|\strong{X*} in a pool with g members is:
#'
#' Y* = g beta_0 + \strong{beta_x}^T \strong{X*} + e*, e* ~ N(0, g sigsq)
#'
#' The assay targets Ybar, the mean Y value for each pool, from which the sum Y*
#' can be calculated as Y* = g Ybar. But the Ybar's may be subject to processing
#' error and/or measurement error. Suppose Ybartilde is the imprecise version of
#' Ybar from the assay. If both errors are present, the assumed error structure
#' is:
#'
#' Ybartilde = Ybar + e_p I(g > 1) + e_m, e_p ~ N(0, sigsq_p),
#' e_m ~ N(0, sigsq_m)
#'
#' with the processing error e_p and measurement error e_m assumed independent
#' of each other. This motivates a maximum likelihood analysis for estimating
#' \strong{theta} = (beta_0, \strong{beta_x}^T)^T based on observed
#' (Ytilde*, \strong{X}*) values, where Ytilde* = g Ytildebar.
#'
#'
#' @param g Numeric vector with pool sizes, i.e. number of members in each pool.
#' @param ytilde Numeric vector (or list of numeric vectors, if some pools have
#' replicates) with poolwise sum \code{Ytilde} values.
#' @param x Numeric matrix with poolwise \strong{\code{X}} values (if any), with
#' one row for each pool. Can be a vector if there is only 1 covariate.
#' @param errors Character string specifying the errors that \code{Y} is subject
#' to. Choices are \code{"neither"}, \code{"processing"} for processing error
#' only, \code{"measurement"} for measurement error only, and \code{"both"}.
#' @param estimate_var Logical value for whether to return variance-covariance
#' matrix for parameter estimates.
#' @param ... Additional arguments to pass to \code{\link[stats]{nlminb}}.
#'
#'
#' @return List of parameter estimates, variance-covariance matrix (if
#' requested), \code{\link[stats]{nlminb}} object, and AIC.
#'
#'
#' @references
#' Schisterman, E.F., Vexler, A., Mumford, S.L. and Perkins, N.J. (2010) "Hybrid
#' pooled-unpooled design for cost-efficient measurement of biomarkers."
#' \emph{Stat. Med.} \strong{29}(5): 597--613.
#'
#'
#' @examples
#' # Load dataset containing data frame with (g, X1*, X2*, Y*, Ytilde*) values
#' # for 500 pools each of size 1, 2, and 3, and list of Ytilde values where 20
#' # of the single-specimen pools have replicates. Ytilde values are affected by
#' # processing error and measurement error; true parameter values are
#' # beta_0 = 0.25, beta_x1 = 0.5, beta_x2 = 0.25, sigsq = 1.
#' data(dat_p_linreg_yerrors)
#' dat <- dat_p_linreg_yerrors$dat
#' reps <- dat_p_linreg_yerrors$reps
#'
#' # Fit Ytilde* vs. (X1*, X2*) ignoring errors in Ytilde (leads to loss of
#' # precision and overestimated sigsq, but no bias).
#' fit.naive <- p_linreg_yerrors(
#'   g = dat$g,
#'   y = dat$y,
#'   x = dat[, c("x1", "x2")],
#'   errors = "neither"
#' )
#' fit.naive$theta.hat
#'
#' # Account for errors in Ytilde*, without using replicates
#' fit.corrected.noreps <- p_linreg_yerrors(
#'   g = dat$g,
#'   y = dat$ytilde,
#'   x = dat[, c("x1", "x2")],
#'   errors = "both"
#' )
#' fit.corrected.noreps$theta.hat
#'
#' # Account for errors in Ytilde*, incorporating the 20 replicates
#' fit.corrected.reps <- p_linreg_yerrors(
#'   g = dat$g,
#'   y = reps,
#'   x = dat[, c("x1", "x2")],
#'   errors = "both"
#' )
#' fit.corrected.reps$theta.hat
#'
#' # In this trial, incorporating replicates resulted in much better estimates
#' # of sigsq (truly 1), sigsq_p (truly 0.4), and sigsq_m (truly = 0.2) but very
#' # similar regression coefficient estimates.
#' fit.corrected.noreps$theta.hat
#' fit.corrected.reps$theta.hat
#'
#'
#' @export
p_linreg_yerrors <- function(g,
                             ytilde,
                             x = NULL,
                             errors = "both",
                             estimate_var = TRUE, ...) {

  # Check that inputs are valid
  if (! errors %in% c("neither", "processing", "measurement", "both")) {
    stop("The input 'errors' should be set to 'neither', 'processing',
         'measurement', or 'both'.")
  }

  # Sample size
  n <- length(ytilde)

  # Get number of X variables (and assign names)
  if (is.null(x)) {
    n.xvars <- 0
    x.varnames <- NULL
  } else {
    x.varname <- deparse(substitute(x))
    if (class(x) != "matrix") {
      x <- as.matrix(x)
    }
    n.xvars <- ncol(x)
    x.varnames <- colnames(x)
    if (is.null(x.varnames)) {
      if (n.xvars == 1) {
        if (length(grep("$", x.varname, fixed = TRUE)) > 0) {
          x.varname <- substr(x.varname,
                              start = which(unlist(strsplit(x.varname, "")) == "$") + 1,
                              stop = nchar(x.varname))
        }
        x.varnames <- x.varname
      } else {
        x.varnames <- paste("x", 1: n.xvars, sep = "")
      }
    }
  }

  # Get number of betas
  n.betas <- 1 + n.xvars

  # Create vector indicating which observations are pools
  Ig <- ifelse(g > 1, 1, 0)

  # Create matrix of (g, X) values
  gx <- cbind(g, x)

  # Separate into replicates and singles
  class.ytilde <- class(ytilde)
  if (class.ytilde == "list") {
    k <- sapply(ytilde, length)
    which.r <- which(k > 1)
    n.r <- length(which.r)
    some.r <- n.r > 0
    if (some.r) {

      # Replicates
      g.r <- g[which.r]
      Ig.r <- Ig[which.r]
      k.r <- k[which.r]
      ytilde.r <- ytilde[which.r]
      gx.r <- gx[which.r, , drop = FALSE]

    }
    n <- n - n.r
    some.s <- n > 0
    if (some.s) {

      # Singles
      g <- g[-which.r]
      Ig <- Ig[-which.r]
      ytilde <- unlist(ytilde[-which.r])
      gx <- gx[-which.r, , drop = FALSE]

    }
  } else {
    some.r <- FALSE
    some.s <- TRUE
  }

  # Get indices for parameters being estimated and create labels
  loc.betas <- 1: n.betas
  beta.labels <- paste("beta", c("0", x.varnames), sep = "_")

  loc.sigsq <- n.betas + 1

  if (errors == "neither") {
    theta.labels <- c(beta.labels, "sigsq")
  } else if (errors == "processing") {
    theta.labels <- c(beta.labels, "sigsq", "sigsq_p")
  } else if (errors == "measurement") {
    theta.labels <- c(beta.labels, "sigsq", "sigsq_m")
  } else if (errors == "both") {
    theta.labels <- c(beta.labels, "sigsq", "sigsq_p", "sigsq_m")
  }

  # Log-likelihood function
  llf <- function(f.theta) {

    # Extract parameters
    f.betas <- matrix(f.theta[loc.betas], ncol = 1)
    f.sigsq <- f.theta[loc.sigsq]

    if (errors == "neither") {
      f.sigsq_p <- 0
      f.sigsq_m <- 0
    } else if (errors == "measurement") {
      f.sigsq_p <- 0
      f.sigsq_m <- f.theta[loc.sigsq + 1]
    } else if (errors == "processing") {
      f.sigsq_p <- f.theta[loc.sigsq + 1]
      f.sigsq_m <- 0
    } else if (errors == "both") {
      f.sigsq_p <- f.theta[loc.sigsq + 1]
      f.sigsq_m <- f.theta[loc.sigsq + 2]
    }

    # Likelihood:
    # L_i = f(Ytilde|X)

    if (some.r) {

      ll.vals <- c()
      for (ii in 1: length(ytilde.r)) {

        # Values for ith subject
        g_i <- g.r[ii]
        Ig_i <- Ig.r[ii]
        k_i <- k.r[ii]
        ytilde_i <- ytilde.r[[ii]]
        gx_i <- gx.r[ii, ]

        # E(Ytilde|X) and V(Ytilde|X)
        Mu_ytilde.x <- matrix(gx_i %*% f.betas, nrow = k_i)
        Sigma_ytilde.x <- g_i * f.sigsq +
          g_i^2 * f.sigsq_p * Ig_i +
          diag(g_i^2 * f.sigsq_m, k_i)

        # Log-likelihood
        ll.vals[ii] <- dmvnorm(x = ytilde_i, log = TRUE,
                               mean = Mu_ytilde.x,
                               sigma = Sigma_ytilde.x)

      }
      ll.r <- sum(ll.vals)

    } else {
      ll.r <- 0
    }

    if (some.s) {

      # E(Ytilde|X) and V(Ytilde|X)
      mu_ytilde.x <- gx %*% f.betas
      sigsq_ytilde.x <- g * f.sigsq + g^2 * f.sigsq_p * Ig + g^2 * f.sigsq_m

      # Log-likelihood
      ll.s <- sum(dnorm(x = ytilde, log = TRUE,
                        mean = mu_ytilde.x, sd = sqrt(sigsq_ytilde.x)))

    } else {
      ll.s <- 0
    }

    # Return negative log-likelihood
    ll <- ll.r + ll.s
    return(-ll)

  }

  # Create list of extra arguments, and assign default starting values and
  # lower values if not specified by user
  extra.args <- list(...)
  if (is.null(extra.args$start)) {
    if (errors == "neither") {
      extra.args$start <- c(rep(0.01, n.betas), 1)
    } else if (errors %in% c("measurement", "processing")) {
      extra.args$start <- c(rep(0.01, n.betas), rep(1, 2))
    } else if (errors == "both") {
      extra.args$start <- c(rep(0.01, n.betas), rep(1, 3))
    }
  }
  if (is.null(extra.args$lower)) {
    if (errors == "neither") {
      extra.args$lower <- c(rep(-Inf, n.betas), 1e-3)
    } else if (errors %in% c("measurement", "processing")) {
      extra.args$lower <- c(rep(-Inf, n.betas), rep(1e-3, 2))
    } else if (errors == "both") {
      extra.args$lower <- c(rep(-Inf, n.betas), rep(1e-3, 3))
    }
  }
  if (is.null(extra.args$control$rel.tol)) {
    extra.args$control$rel.tol <- 1e-6
  }
  if (is.null(extra.args$control$eval.max)) {
    extra.args$control$eval.max <- 1000
  }
  if (is.null(extra.args$control$iter.max)) {
    extra.args$control$iter.max <- 750
  }

  # Obtain ML estimates
  ml.max <- do.call(nlminb, c(list(objective = llf), extra.args))

  # Create list to return
  theta.hat <- ml.max$par
  names(theta.hat) <- theta.labels
  ret.list <- list(theta.hat = theta.hat)

  # If requested, add variance-covariance matrix to ret.list
  if (estimate_var) {

    hessian.mat <- hessian(f = llf, x0 = theta.hat)
    theta.variance <- try(solve(hessian.mat), silent = TRUE)
    if (class(theta.variance) == "try-error") {
      message("Estimated Hessian matrix is singular, so variance-covariance matrix cannot be obtained.")
      ret.list$theta.var <- NULL
    } else {
      colnames(theta.variance) <- rownames(theta.variance) <- theta.labels
      ret.list$theta.var <- theta.variance
    }

  }

  # Add nlminb object and AIC to ret.list
  ret.list$nlminb.object <- ml.max
  ret.list$aic <- 2 * (length(theta.hat) + ml.max$objective)

  # Return ret.list
  return(ret.list)

}
