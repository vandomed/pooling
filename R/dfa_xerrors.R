#' Discriminant Function Approach for Estimating Odds Ratio with Normal Exposure
#' Subject to Measurement Error
#'
#' Assumes exposure measurements are subject to additive normal measurement
#' error, and exposure given covariates and outcome is a normal-errors linear
#' regression. Some replicates are required for identifiability.
#'
#'
#' @inheritParams logreg_xerrors
#' @param merror Logical value for whether there is measurement error.
#' @param ... Additional arguments to pass to \code{\link[stats]{nlminb}}.
#'
#'
#' @return
#' List containing:
#' \enumerate{
#'  \item Numeric vector with point estimates for \code{gamma_y}, \code{sigsq},
#'  and the covariate-adjusted log-odds ratio, and the estimated variance for
#'  the log-odds ratio estimate if \code{estimate_var = TRUE}.
#'  \item Object returned by \code{\link[stats]{nlminb}} containing information
#'  about likelihood maximization.
#'  \item Akaike information criterion (AIC).
#' }
#'
#'
#' @references
#' Lyles, R.H., Van Domelen, D.R., Mitchell, E.M. and Schisterman, E.F. (2015)
#' "A discriminant function approach to adjust for processing and measurement
#' error When a biomarker is assayed in pooled samples."
#' \emph{Int. J. Environ. Res. Public Health} \strong{12}(11): 14723--14740.
#'
#'
#' @examples
#' # Load dataset - dat1 has (Y, C) values and dat1.xtilde is list with 1 or 2
#' # Xtilde measurements for each subject.
#' data(dat1)
#' data(dat1_xtilde)
#'
#' # Estimate log-OR for X and Y adjusted for C, ignoring measurement error
#' fit1 <- dfa_xerrors(y = dat1$y, xtilde = dat1_xtilde, c = dat1$c, merror = FALSE)
#' fit1$estimates
#'
#' # Repeat, but accounting for measurement error. Closer to true log-OR of 0.5.
#' fit2 <- dfa_xerrors(y = dat1$y, xtilde = dat1_xtilde, c = dat1$c, merror = TRUE)
#' fit2$estimates
#'
#'
#' @export
dfa_xerrors <- function(y, xtilde, c = NULL, merror = TRUE, ...) {

  # Get name of y input
  y.varname <- deparse(substitute(y))
  if (length(grep("$", y.varname, fixed = TRUE)) > 0) {
    y.varname <- substr(y.varname,
                        start = which(unlist(strsplit(y.varname, "")) == "$") + 1,
                        stop = nchar(y.varname))
  }

  # Get number of C variables (and assign names)
  if (is.null(c)) {
    c.varnames <- NULL
    n.cvars <- 0
  } else {
    c.varname <- deparse(substitute(c))
    if (class(c) != "matrix") {
      c <- as.matrix(c)
    }
    n.cvars <- ncol(c)
    c.varnames <- colnames(c)
    if (is.null(c.varnames)) {
      if (n.cvars == 1) {
        if (length(grep("$", c.varname, fixed = TRUE)) > 0) {
          c.varname <- substr(c.varname,
                              start = which(unlist(strsplit(c.varname, "")) == "$") + 1,
                              stop = nchar(c.varname))
        }
        c.varnames <- c.varname
      } else {
        c.varnames <- paste("c", 1: n.cvars, sep = "")
      }
    }
  }

  # Get number of gammas
  n.gammas <- 2 + n.cvars

  # Sample size
  n <- length(y)

  # Construct (1, Y, C) matrix
  oneyc <- cbind(rep(1, n), y, c)

  # If no measurement error and xtilde is a list, just use first measurements
  if (! merror & class(xtilde) == "list") {
    xtilde <- sapply(xtilde, function(x) x[1])
  }

  # Separate out replicates
  class.xtilde <- class(xtilde)
  if (class.xtilde == "list") {
    k <- sapply(xtilde, length)
    which.r <- which(k > 1)
    n.r <- length(which.r)
    some.r <- n.r > 0
    if (some.r) {

      # Replicates
      k.r <- k[which.r]
      y.r <- y[which.r]
      oneyc.r <- oneyc[which.r, , drop = FALSE]
      xtilde.r <- xtilde[which.r]

    }
    n <- n - n.r
    some.s <- n > 0
    if (some.s) {

      # Singles
      y <- y[-which.r]
      oneyc <- oneyc[-which.r, , drop = FALSE]
      xtilde <- unlist(xtilde[-which.r])

    }
  } else {
    some.r <- FALSE
    some.s <- TRUE
  }

  # Get indices for parameters being estimated and create labels
  loc.gammas <- 1: n.gammas
  gamma.labels <- paste("gamma", c("0", y.varname, c.varnames), sep = "_")

  loc.sigsq <- n.gammas + 1
  loc.sigsq_m <- n.gammas + 2

  theta.labels <- c(gamma.labels, "sigsq")
  if (merror) {
    theta.labels <- c(theta.labels, "sigsq_m")
  }

  # Log-likelihood function
  ll.f <- function(f.theta) {

    # Extract parameters
    f.gammas <- matrix(f.theta[loc.gammas], ncol = 1)
    f.sigsq <- f.theta[loc.sigsq]

    if (merror) {

      # Likelihood:
      # L_i = f(Xtilde|Y,C)

      f.sigsq_m <- f.theta[loc.sigsq_m]

      if (some.r) {

        ll.vals <- c()
        for (ii in 1: length(xtilde.r)) {

          # Values for ith subject
          k_i <- k.r[ii]
          oneyc_i <- oneyc.r[ii, ]
          xtilde_i <- xtilde.r[[ii]]

          # E(Xtilde|Y,C) and V(Xtilde|Y,C)
          Mu_xtilde.yc <- matrix(oneyc_i %*% f.gammas, ncol = k_i)
          Sigma_xtilde.yc <- f.sigsq + diag(x = f.sigsq_m, ncol = k_i, nrow = k_i)

          # Log-likelihood
          ll.vals[ii] <- dmvnorm(x = xtilde_i, log = TRUE,
                                 mean = Mu_xtilde.yc,
                                 sigma = Sigma_xtilde.yc)

        }
        ll.r <- sum(ll.vals)

      } else {
        ll.r <- 0
      }

      if (some.s) {

        # E(Xtilde|Y,C)
        mu_xtilde.yc <- oneyc %*% f.gammas

        # Log-likelihood
        ll.s <- sum(dnorm(x = xtilde, log = TRUE,
                          mean = mu_xtilde.yc,
                          sd = sqrt(f.sigsq + f.sigsq_m)))

      } else {
        ll.s <- 0
      }

      ll <- ll.r + ll.s

    } else {

      # Likelihood:
      # L = f(X|Y,C)

      # E(X|Y,C)
      mu_x.yc <- oneyc %*% f.gammas

      # Log-likelihood
      ll <- sum(dnorm(x = xtilde, log = TRUE, mean = mu_x.yc, sd = sqrt(f.sigsq)))

    }

    # Return negative log-likelihood
    return(-ll)

  }

  # Create list of extra arguments, and assign default starting values and
  # lower values if not specified by user
  extra.args <- list(...)
  if (is.null(extra.args$start)) {
    if (! merror) {
      extra.args$start <- c(rep(0.01, n.gammas), 1)
    } else {
      extra.args$start <- c(rep(0.01, n.gammas), 1, 1)
    }
  }
  if (is.null(extra.args$lower)) {
    if (! merror) {
      extra.args$lower <- c(rep(-Inf, n.gammas), 1e-3)
    } else {
      extra.args$lower <- c(rep(-Inf, n.gammas), rep(1e-3, 2))
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
  ml.max <- do.call(nlminb, c(list(objective = ll.f), extra.args))
  ml.estimates <- ml.max$par

  # Obtain point estimate for log-odds ratio
  gamma_y.hat <- ml.estimates[2]
  sigsq.hat <- ml.estimates[loc.sigsq]
  logOR.hat <- gamma_y.hat / sigsq.hat

  # Estimate variance of logOR.hat and perform bias adjustment
  hessian.mat <- pracma::hessian(f = ll.f, x0 = ml.estimates)
  theta.variance <- try(solve(hessian.mat), silent = TRUE)
  if (class(theta.variance) == "try-error") {
    message("Estimated Hessian matrix is singular, so variance-covariance matrix cannot be obtained and bias adjustment cannot be applied.")
    theta.variance <- NULL
    logOR.var <- logOR.var <- logOR_adj.hat <- logOR_adj.var <- NA
  } else {
    fprime <- matrix(c(1 / sigsq.hat, -gamma_y.hat / sigsq.hat^2), nrow = 1)
    colnames(theta.variance) <- rownames(theta.variance) <- theta.labels
    logOR.var <- fprime %*%
      theta.variance[c(2, loc.sigsq), c(2, loc.sigsq)] %*% t(fprime)
    sigsq.var <- theta.variance[loc.sigsq, loc.sigsq]
    logOR_adj.hat <- logOR.hat - gamma_y.hat * sigsq.var / sigsq.hat^3
    logOR_adj.var <- logOR.var * (logOR_adj.hat / logOR.hat)^2
    if (sign(logOR.hat) != sign(logOR_adj.hat)) {
      message("Bias adjustment flipped the sign of the log-OR estimate, so you may want to use the non-bias adjusted version.")
    }
  }

  # Create vector of estimates to return
  estimates <- c(ml.estimates,
                 logOR.hat, logOR.var,
                 logOR_adj.hat, logOR_adj.var)
  names(estimates) <- c(theta.labels,
                        "logOR.hat", "logOR.var",
                        "logOR_adj.hat", "logOR_adj.var")

  # Create list of estimates, returned nlminb object, and AIC to return
  ret.list <- list(estimates = estimates,
                   theta.var = theta.variance,
                   nlminb.object = ml.max,
                   aic = 2 * (length(ml.estimates) + ml.max$objective))
  return(ret.list)

}
