#' Discriminant Function Approach for Estimating Odds Ratio with Normal Exposure
#' Measured in Pools and Subject to Errors
#'
#' Assumes exposure measurements are subject to additive normal processing error
#' and measurement error, and exposure given covariates and outcome is a
#' normal-errors linear regression.
#'
#'
#' @inheritParams p_logreg_xerrors
#'
#' @param y Numeric vector of poolwise \code{Y} values (number of cases in each
#' pool).
#'
#' @param constant_or Logical value for whether to assume a constant OR for
#' \code{X}, which means that \code{sigsq_1 = sigsq_0}. If \code{NULL}, model is
#' fit with and without this assumption.
#'
#'
#' @return
#' List of point estimates, variance-covariance matrix, objects returned by
#' \code{\link[stats]{nlminb}}, and AICs, for one or two models depending on
#' \code{constant_or}.
#'
#'
#' @references
#' Lyles, R.H., Van Domelen, D.R., Mitchell, E.M. and Schisterman, E.F. (2015)
#' "A discriminant function approach to adjust for processing and measurement
#' error When a biomarker is assayed in pooled samples."
#' \emph{Int. J. Environ. Res. Public Health} \strong{12}(11): 14723--14740.
#'
#' Schisterman, E.F., Vexler, A., Mumford, S.L. and Perkins, N.J. (2010) "Hybrid
#' pooled-unpooled design for cost-efficient measurement of biomarkers."
#' \emph{Stat. Med.} \strong{29}(5): 597--613.
#'
#'
#' @export
p_dfa_xerrors <- function(g, y, xtilde, c = NULL,
                          constant_or = NULL,
                          errors = "both", ...) {

  # Check that inputs are valid
  if (! is.null(constant_or) && ! is.logical(constant_or)) {
    stop("The input 'contant_or' should be set to TRUE, FALSE, or NULL.")
  }
  if (! errors %in% c("neither", "processing", "measurement", "both")) {
    stop("The input 'errors' should be set to 'neither', 'processing',
         'measurement', or 'both'.")
  }

  # Sample size
  n <- length(y)

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

  # Separate into replicate and singles
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
      Ig.r <- Ig[which.r]
      y.r <- y[which.r]
      gyc.r <- gyc[which.r, , drop = FALSE]
      xtilde.r <- xtilde[which.r]

    }
    n <- n - n.r
    some.s <- n > 0
    if (some.s) {

      # Singles
      g.s <- g[-which.r]
      Ig.s <- Ig[-which.r]
      y.s <- y[-which.r]
      gyc.s <- gyc[-which.r, , drop = FALSE]
      xtilde.s <- unlist(xtilde[-which.r])

    }
  } else {
    some.r <- FALSE
    some.s <- TRUE
  }

  # Get indices for parameters being estimated and create labels
  loc.gammas <- 1: n.gammas
  gamma.labels <- paste("gamma", c("0", y.varname, c.varnames), sep = "_")

  loc.sigsq <- n.gammas + 1
  loc.sigsq_1 <- n.gammas + 1
  loc.sigsq_0 <- n.gammas + 2

  if (errors == "neither") {
    theta.labels1 <- c(gamma.labels, "sigsq_1", "sigsq_0")
    theta.labels2 <- c(gamma.labels, "sigsq")
  } else if (errors == "processing") {
    theta.labels1 <- c(gamma.labels, "sigsq_1", "sigsq_0", "sigsq_p")
    theta.labels2 <- c(gamma.labels, "sigsq", "sigsq_p")
  } else if (errors == "measurement") {
    theta.labels1 <- c(gamma.labels, "sigsq_1", "sigsq_0", "sigsq_m")
    theta.labels2 <- c(gamma.labels, "sigsq", "sigsq_m")
  } else if (errors == "both") {
    theta.labels1 <- c(gamma.labels, "sigsq_1", "sigsq_0", "sigsq_p", "sigsq_m")
    theta.labels2 <- c(gamma.labels, "sigsq", "sigsq_p", "sigsq_m")
  }

  # Fit model with different residual error variances
  if (is.null(constant_or) || ! constant_or) {

    # Log-likelihood function
    ll.f1 <- function(f.theta) {

      # Extract parameters
      f.gammas <- matrix(f.theta[loc.gammas], ncol = 1)
      f.sigsq_1 <- f.theta[loc.sigsq_1]
      f.sigsq_0 <- f.theta[loc.sigsq_0]

      if (errors == "neither") {
        f.sigsq_p <- 0
        f.sigsq_m <- 0
      } else if (errors == "measurement") {
        f.sigsq_p <- 0
        f.sigsq_m <- f.theta[loc.sigsq_0 + 1]
      } else if (errors == "processing") {
        f.sigsq_p <- f.theta[loc.sigsq_0 + 1]
        f.sigsq_m <- 0
      } else if (errors == "both") {
        f.sigsq_p <- f.theta[loc.sigsq_0 + 1]
        f.sigsq_m <- f.theta[loc.sigsq_0 + 2]
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
          y_i <- y.r[ii]
          gyc_i <- gyc.r[ii, ]
          xtilde_i <- xtilde.r[[ii]]

          # E(Xtilde|Y,C) and V(Xtilde|Y,C)
          Mu_xtilde.yc <- matrix(gyc_i %*% f.gammas, ncol = k_i)
          Sigma_xtilde.yc <-
            matrix(g_i * ifelse(y_i, f.sigsq_1, f.sigsq_0) +
                     g_i^2 * f.sigsq_p * Ig_i, ncol = k_i, nrow = k_i) +
            diag(x = g_i^2 * f.sigsq_m, ncol = k_i, nrow = k_i)

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

        # E(Xtilde|Y,C) and V(Xtilde|Y,C)
        mu_xtilde.yc <- gyc.s %*% f.gammas
        sigsq_xtilde.yc <- g.s * ifelse(y.s, f.sigsq_1, f.sigsq_0) +
          g.s^2 * f.sigsq_p * Ig.s + g.s^2 * f.sigsq_m

        # Log-likelihood
        ll.s <- sum(dnorm(x = xtilde.s, log = TRUE,
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
        extra.args$start <- c(rep(0.01, n.gammas), rep(1, 2))
      } else if (errors %in% c("measurement", "processing")) {
        extra.args$start <- c(rep(0.01, n.gammas), rep(1, 3))
      } else if (errors == "both") {
        extra.args$start <- c(rep(0.01, n.gammas), rep(1, 4))
      }
    }
    if (is.null(extra.args$lower)) {
      if (errors == "neither") {
        extra.args$lower <- c(rep(-Inf, n.gammas), rep(1e-3, 2))
      } else if (errors %in% c("measurement", "processing")) {
        extra.args$lower <- c(rep(-Inf, n.gammas), rep(1e-3, 3))
      } else if (errors == "both") {
        extra.args$lower <- c(rep(-Inf, n.gammas), rep(1e-3, 4))
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
    ml.max1 <- do.call(nlminb, c(list(objective = ll.f1), extra.args))
    ml.estimates1 <- ml.max1$par

    # Variance estimates
    hessian.mat1 <- pracma::hessian(f = ll.f1, x0 = ml.estimates1)
    theta.variance1 <- try(solve(hessian.mat1), silent = TRUE)
    if (class(theta.variance1) == "try-error") {
      message("Estimated Hessian matrix is singular, so variance-covariance matrix cannot be obtained.")
      theta.variance1 <- NULL
    } else {
      colnames(theta.variance1) <- rownames(theta.variance1) <- theta.labels1
    }

    # Create vector of estimates and calculate AIC
    estimates1 <- ml.estimates1
    names(estimates1) <- theta.labels1
    theta.var1 <- theta.variance1
    aic1 <- 2 * (length(estimates1) + ml.max1$objective)

  }

  # Fit model with same residual error variances
  if (is.null(constant_or) || constant_or) {

    # Log-likelihood function
    ll.f2 <- function(f.theta) {

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
          y_i <- y.r[ii]
          gyc_i <- gyc.r[ii, ]
          xtilde_i <- xtilde.r[[ii]]

          # E(Xtilde|Y,C) and V(Xtilde|Y,C)
          Mu_xtilde.yc <- matrix(gyc_i %*% f.gammas, ncol = k_i)
          Sigma_xtilde.yc <-
            matrix(g_i * f.sigsq +
                     g_i^2 * f.sigsq_p * Ig_i, ncol = k_i, nrow = k_i) +
            diag(x = g_i^2 * f.sigsq_m, ncol = k_i, nrow = k_i)

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

        # E(Xtilde|Y,C) and V(Xtilde|Y,C)
        mu_xtilde.yc <- gyc.s %*% f.gammas
        sigsq_xtilde.yc <- g.s * f.sigsq +
          g.s^2 * f.sigsq_p * Ig.s +
          g.s^2 * f.sigsq_m

        # Log-likelihood
        ll.s <- sum(dnorm(x = xtilde.s, log = TRUE,
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
        extra.args$start <- c(rep(0.01, n.gammas), rep(1, 2))
      } else if (errors == "both") {
        extra.args$start <- c(rep(0.01, n.gammas), rep(1, 3))
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
    ml.max2 <- do.call(nlminb, c(list(objective = ll.f2), extra.args))
    ml.estimates <- ml.max2$par

    # Obtain point estimate for log-odds ratio
    gamma_y.hat <- ml.estimates[2]
    sigsq.hat <- ml.estimates[loc.sigsq]
    logOR.hat <- gamma_y.hat / sigsq.hat

    # Estimate variance of logOR.hat and perform bias adjustment
    hessian.mat <- pracma::hessian(f = ll.f2, x0 = ml.estimates)
    theta.variance <- try(solve(hessian.mat), silent = TRUE)
    if (class(theta.variance) == "try-error") {
      message("Estimated Hessian matrix is singular, so variance-covariance matrix cannot be obtained and bias adjustment cannot be applied.")
      theta.variance <- NULL
      logOR.var <- logOR.var <- logOR_adj.hat <- logOR_adj.var <- NA
    } else {
      fprime <- matrix(c(1 / sigsq.hat, -gamma_y.hat / sigsq.hat^2), nrow = 1)
      colnames(theta.variance) <- rownames(theta.variance) <- theta.labels2
      logOR.var <- fprime %*%
        theta.variance[c(2, loc.sigsq), c(2, loc.sigsq)] %*% t(fprime)
      sigsq.var <- theta.variance[loc.sigsq, loc.sigsq]
      logOR_adj.hat <- logOR.hat - gamma_y.hat * sigsq.var / sigsq.hat^3
      logOR_adj.var <- logOR.var * (logOR_adj.hat / logOR.hat)^2
      if (sign(logOR.hat) != sign(logOR_adj.hat)) {
        message("Bias adjustment flipped the sign of the log-OR estimate, so you may want to use the non-bias adjusted version.")
      }
    }

    # Create vector of estimates and calculate AIC
    estimates2 <- c(ml.estimates, logOR.hat, logOR.var, logOR_adj.hat, logOR_adj.var)
    names(estimates2) <- c(theta.labels2, "logOR.hat", "logOR.var", "logOR_adj.hat", "logOR_adj.var")
    theta.var2 <- theta.variance
    aic2 <- 2 * (length(ml.estimates) + ml.max2$objective)

  }

  # Return objects
  if (is.null(constant_or)) {
    return(list(estimates1 = estimates1,
                estimates2 = estimates2,
                theta.var1 = theta.var1,
                theta.var2 = theta.var2,
                nlminb.object1 = ml.max1,
                nlminb.object2 = ml.max2,
                aic1 = aic1,
                aic2 = aic2))
  } else if (constant_or) {
    return(list(estimates = estimates2,
                theta.var = theta.var2,
                nlminb.object = ml.max2,
                aic = aic2))
  } else {
    return(list(estimates = estimates1,
                theta.var = theta.var1,
                nlminb.object = ml.max1,
                aic = aic1))
  }
}
