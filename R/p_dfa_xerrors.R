#' Discriminant Function Approach for Estimating Odds Ratio with Poolwise Data
#' and Errors in Continuous Exposure
#'
#' Implements method described in [1].
#'
#'
#' @inheritParams p_logreg_xerrors
#' @param y Numeric vector of poolwise \code{Y} values (number of cases in each
#' pool).
#'
#'
#' @return List containing:
#' \enumerate{
#'   \item Named numeric vector with point estimates for \code{gamma_y},
#'   \code{sigsq}, and the covariate-adjusted log-odds ratio, and the estimated
#'  variance for the log-odds ratio estimate if \code{estimate.var = TRUE}.
#'   \item Object returned by \code{\link[stats]{nlminb}} containing information
#'   about likelihood maximization.
#'   \item Akaike information criterion (AIC).
#'  }
#'
#'
#' @references
#' Lyles, R.H., Van Domelen, D.R., Mitchell, E.M. and Schisterman, E.F. (2015)
#' "A Discriminant Function Approach to Adjust for Processing and Measurement
#' Error When a Biomarker is Assayed in Pooled Samples."
#' \emph{Int. J. Environ. Res. Public Health} \strong{12}(11): 14723--14740.
#'
#' Acknowledgment: This material is based upon work supported by the National
#' Science Foundation Graduate Research Fellowship under Grant No. DGE-0940903.
#'
#'
#' @export
p_dfa_xerrors <- function(g, y, xtilde, c = NULL,
                          errors = "both", ...) {

  # Check that inputs are valid
  if (! errors %in% c("neither", "processing", "measurement", "both")) {
    stop("The input 'errors' should be set to 'neither', 'processing',
         'measurement', or 'both'.")
  }

  # Sample size
  n.pools <- length(y)

  # Get name of y input
  y.varname <- deparse(substitute(y))

  # Get number of C variables (and assign names)
  if (is.null(c)) {
    c.varnames <- NULL
    n.cvars <- 0
  } else {
    if (class(c) != "matrix") {
      c <- as.matrix(c)
    }
    n.cvars <- ncol(c)
    c.varnames <- colnames(c)
    if (is.null(c.varnames)) {
      if (n.cvars == 1) {
        c.varnames <- deparse(substitute(c))
      } else {
        c.varnames <- paste("c", 1: n.cvars, sep = "")
      }
    }
  }

  # Get number of gammas
  n.gammas <- 2 + n.cvars

  # Create vector indicating which observations are pools
  Ig <- ifelse(g > 1, 1, 0)

  # Create matrix of (g, Y, C) values
  gyc <- cbind(g, y, c)

  # Separate out replicate from singles
  class.xtilde <- class(xtilde)
  if (class.xtilde == "list") {
    k <- sapply(xtilde, length)
    which.r <- which(k > 1)
    n.r <- length(which.r)
    some.r <- n.r > 0
    if (some.r) {

      # Replicates
      k.r <- k[which.r]
      g.r <- g[which.r]
      y.r <- y[which.r]
      gyc.r <- gyc[which.r, , drop = FALSE]
      xtilde.r <- xtilde[which.r]

    }
    n <- n - n.r
    some.s <- n > 0
    if (some.s) {

      # Singles
      g <- g[-which.r]
      Ig <- Ig[-which.r]
      y <- y[-which.r]
      gyc <- gyc[-which.r, , drop = FALSE]
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

  if (errors == "neither") {
    theta.labels <- c(gamma.labels, "sigsq")
  } else if (errors == "processing") {
    theta.labels <- c(gamma.labels, "sigsq", "sigsq_p")
  } else if (errors == "measurement") {
    theta.labels <- c(gamma.labels, "sigsq", "sigsq_m")
  } else if (errors == "both") {
    theta.labels <- c(gamma.labels, "sigsq", "sigsq_p", "sigsq_m")
  }

  # Log-likelihood function
  ll.f <- function(f.theta) {

    # Extract parameters
    f.gammas <- matrix(f.theta[loc.gammas], ncol = 1)
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
    # L = f(Xtilde|Y,C)

    if (some.r) {

      ll.vals <- c()
      for (ii in 1: length(xtilde.r)) {

        # Values for ith subject
        g_i <- g.r[ii]
        Ig_i <- Ig.r[ii]
        k_i <- k.r[ii]
        gyc_i <- gyc.r[ii, ]
        xtilde_i <- xtilde.r[[ii]]

        # E(Xtilde|Y,C) and V(Xtilde|Y,C)
        Mu_xtilde.yc <- matrix(gyc_i %*% f.gammas, ncol = k_i)
        Sigma_xtilde.yc <-
          matrix(g_i * f.sigsq + g_i^2 * f.sigsq_p * Ig_i,
                 ncol = k_i, nrow = k_i) +
          diag(x = g_i^2 * f.sigsq_m, ncol = k_i, nrow = k_i)

        # Log-likelihood
        ll.vals[ii] <- mvtnorm::dmvnorm(x = xtilde_i, log = TRUE,
                                        mean = Mu_xtilde.yc,
                                        sigma = Sigma_xtilde.yc)

      }
      ll.r <- sum(ll.vals)

    } else {
      ll.r <- 0
    }

    if (some.s) {

      # E(Xtilde|Y,C) and V(Xtilde|Y,C)
      mu_xtilde.yc <- gyc %*% f.gammas
      sigsq_xtilde.yc <- g * f.sigsq + g^2 * f.sigsq_p * Ig + g^2 * f.sigsq_m

      # Log-likelihood
      ll.s <- sum(dnorm(x = xtilde, log = TRUE,
                        mean = mu_xtilde.yc, sd = sqrt(sigsq_xtilde.yc)))

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
      extra.args$start <- c(rep(0.01, n.gammas), 1)
    } else if (errors %in% c("measurement", "processing")) {
      extra.args$start <- c(rep(0.01, n.gammas), 1, 1)
    } else if (errors == "both") {
      extra.args$start <- c(rep(0.01, n.gammas), 1, 1, 1)
    }
  }
  if (is.null(extra.args$lower)) {
    if (errors == "neither") {
      extra.args$lower <- c(rep(-Inf, n.gammas), 1e-3)
    } else if (errors %in% c("measurement", "processing")) {
      extra.args$lower <- c(rep(-Inf, n.gammas), rep(1e-3, 2))
    } else if (errors == "both") {
      extra.args$lower <- c(rep(-Inf, n.gammas), rep(1e-3, 3))
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
