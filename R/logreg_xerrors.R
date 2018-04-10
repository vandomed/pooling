#' Logistic Regression with Normal Exposure Subject to Errors
#'
#' Assumes exposure measurements are subject to additive normal measurement
#' error, and exposure given covariates is a normal-errors linear regression.
#' Some replicates are required for identifiability.
#'
#'
#' @param xtilde List of numeric vectors with \code{Xtilde} values.
#'
#' @param c Numeric matrix with \strong{\code{C}} values (if any), with
#' one row for each pool. Can be a vector if there is only 1 covariate.
#'
#' @param prev Numeric value specifying disease prevalence, allowing for valid
#' estimation of the intercept with case-control sampling. Can specify
#' \code{samp_y1y0} instead if sampling rates are known.
#'
#' @param samp_y1y0 Numeric vector of length 2 specifying sampling probabilities
#' for cases and controls, allowing for valid estimation of the intercept with
#' case-control sampling. Can specify \code{prev} instead if it's easier.
#'
#' @param merror Logical value for whether there is measurement error.
#'
#' @param approx_integral Logical value for whether to use the probit
#' approximation for the logistic-normal integral, to avoid numerically
#' integrating \code{X}'s out of the likelihood function.
#'
#' @param integrate_tol Numeric value specifying the \code{tol} input to
#' \code{\link{adaptIntegrate}}. Only used if \code{approx_integral = FALSE}.
#'
#' @param integrate_tol_start Same as \code{integrate_tol}, but applies only to
#' the very first iteration of ML maximization. The first iteration tends to
#' take much longer than subsequent ones, so less precise integration at the
#' start can speed things up.
#'
#' @param integrate_tol_hessian Same as \code{integrate_tol}, but for use when
#' estimating the Hessian matrix only. Sometimes more precise integration
#' (i.e. smaller tolerance) than used for maximizing the likelihood helps
#' prevent cases where the inverse Hessian is not positive definite.
#'
#' @param estimate_var Logical value for whether to return variance-covariance
#' matrix for parameter estimates.
#'
#' @param ... Additional arguments to pass to \code{\link[stats]{nlminb}}.
#'
#'
#' @return List containing:
#' \enumerate{
#' \item Numeric vector of parameter estimates.
#' \item Variance-covariance matrix (if \code{estimate_var = TRUE}).
#' \item Returned \code{\link[stats]{nlminb}} object from maximizing the
#' log-likelihood function.
#' \item Akaike information criterion (AIC).
#' }
#'
#'
#' @export
logreg_xerrors <- function(y, xtilde, c = NULL,
                           prev = NULL, samp_y1y0 = NULL,
                           merror = TRUE,
                           approx_integral = TRUE,
                           integrate_tol = 1e-8,
                           integrate_tol_start = integrate_tol,
                           integrate_tol_hessian = integrate_tol,
                           estimate_var = FALSE,
                           ...) {

  # Get name of xtilde input
  x.varname <- deparse(substitute(xtilde))

  # Get information about covariates C
  if (is.null(c)) {
    c.varnames <- NULL
    n.cvars <- 0
    some.cs <- FALSE
  } else {
    c.varname <- deparse(substitute(c))
    if (class(c) != "matrix") {
      c <- as.matrix(c)
    }
    n.cvars <- ncol(c)
    some.cs <- TRUE
    c.varnames <- colnames(c)
    if (is.null(c.varnames)) {
      if (n.cvars == 1) {
        c.varnames <- c.varname
      } else {
        c.varnames <- paste("c", 1: n.cvars, sep = "")
      }
    }
  }

  # Get number of betas and alphas
  n.betas <- 2 + n.cvars
  n.alphas <- 1 + n.cvars

  # Sample sizes
  n <- length(y)
  locs.cases <- which(y == 1)
  locs.controls <- which(y == 0)
  n1 <- length(locs.cases)
  n0 <- length(locs.controls)

  # Construct (1, C) matrix
  onec <- cbind(rep(1, n), c)

  # Calculate offsets if prev or samp_y1y0 specified
  if (! is.null(prev)) {
    q <- rep(log(n1 / n0 * prev / (1 - prev)), n)
  } else if (! is.null(samp_y1y0)) {
    q <- rep(log(samp_y1y0[1] / samp_y1y0[2]), n)
  } else {
    q <- rep(0, n)
  }

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
      y.r <- y[which.r]
      c.r <- c[which.r, , drop = FALSE]
      q.r <- q[which.r]
      onec.r <- onec[which.r, , drop = FALSE]
      xtilde.r <- xtilde[which.r]

    }
    n <- n - n.r
    some.s <- n > 0
    if (some.s) {

      # Singles
      k <- k[-which.r]
      y <- y[-which.r]
      c <- c[-which.r, , drop = FALSE]
      q <- q[-which.r]
      onec <- onec[-which.r, , drop = FALSE]
      xtilde <- unlist(xtilde[-which.r])

    }
  } else {
    some.r <- FALSE
    some.s <- TRUE
  }

  # Create (1, X, C) matrix if merror is FALSE
  if (! merror) {
    onexc <- cbind(rep(1, n), xtilde, c)
  }

  # Get indices for parameters being estimated and create labels
  loc.betas <- 1: n.betas
  beta.labels <- paste("beta", c("0", x.varname, c.varnames), sep = "_")

  loc.alphas <- (n.betas + 1): (n.betas + n.alphas)
  alpha.labels <- paste("alpha", c("0", c.varnames), sep = "_")

  loc.sigsq_x.c <- n.betas + n.alphas + 1
  loc.sigsq_m <- n.betas + n.alphas + 2

  theta.labels <- c(beta.labels, alpha.labels, "sigsq_x.c")
  if (merror) {
    theta.labels <- c(theta.labels, "sigsq_m")
  }

  # Log-likelihood function
  ll.f <- function(f.theta, estimating.hessian = FALSE) {

    # Extract parameters
    f.betas <- matrix(f.theta[loc.betas], ncol = 1)
    f.beta_0 <- f.betas[1]
    f.beta_x <- f.betas[2]
    f.beta_c <- matrix(f.betas[-c(1: 2)], ncol = 1)

    f.alphas <- matrix(f.theta[loc.alphas], ncol = 1)
    f.alpha_0 <- f.alphas[1]
    f.alpha_c <- matrix(f.alphas[-1], ncol = 1)

    f.sigsq_x.c <- f.theta[loc.sigsq_x.c]

    if (merror) {

      # Likelihood (Xtilde vector-valued for replicates)
      # L = f(Y, Xtilde|C)
      #   = [\int_X f(Y|X,C) f(X|Xtilde,C) dX] f(Xtilde|C)
      #   = int_X f(Y|X,C) f(Xtilde|X) f(X|C) dX

      f.sigsq_m <- f.theta[loc.sigsq_m]

      # Set skip.rest to FALSE (want to break from loop if get 0 integral)
      skip.rest <- FALSE

      if (some.r) {

        # E(X|C)
        mu_x.c <- onec.r %*% f.alphas

        if (approx_integral) {

          # Probit approximation for logistic-normal integral

          ll.vals <- c()
          for (ii in 1: length(xtilde.r)) {

            # Values for ith subject
            k_i <- k.r[ii]
            y_i <- y.r[ii]
            c_i <- c.r[ii, ]
            q_i <- q.r[ii]
            mu_x.c_i <- mu_x.c[ii]
            xtilde_i <- xtilde.r[[ii]]

            # E(X|Xtilde,C) and V(X|Xtilde,C)
            Mu_xxtilde.c <- matrix(mu_x.c_i, nrow = k_i + 1)
            Sigma_xxtilde.c_11 <- f.sigsq_x.c
            Sigma_xxtilde.c_12 <- matrix(f.sigsq_x.c, ncol = k_i)
            Sigma_xxtilde.c_21 <- t(Sigma_xxtilde.c_12)
            Sigma_xxtilde.c_22 <- f.sigsq_x.c + diag(x = f.sigsq_m, ncol = k_i, nrow = k_i)

            mu_x.xtildec <- Mu_xxtilde.c[1] + Sigma_xxtilde.c_12 %*%
              solve(Sigma_xxtilde.c_22) %*% (xtilde_i - Mu_xxtilde.c[-1])
            sigsq_x.xtildec <- Sigma_xxtilde.c_11 - Sigma_xxtilde.c_12 %*%
              solve(Sigma_xxtilde.c_22) %*% Sigma_xxtilde.c_21

            # Approximation of \int_X f(Y|X,C) f(X|Xtilde,C) dX
            if (some.cs) {
              t <- (f.beta_0 + f.beta_x * mu_x.xtildec +
                      c_i %*% f.beta_c + q_i) /
                sqrt(1 + sigsq_x.xtildec * f.beta_x^2 / 1.7^2)
            } else {
              t <- (f.beta_0 + f.beta_x * mu_x.xtildec + q_i) /
                sqrt(1 + sigsq_x.xtildec * f.beta_x^2 / 1.7^2)
            }
            p <- exp(t) / (1 + exp(t))
            part1 <- dbinom(x = y_i, size = 1, prob = p, log = TRUE)

            # log[f(Xtilde|C)]
            if (k_i == 2) {

              mu_xtilde1.xtilde2c <- mu_x.c_i + Sigma_xxtilde.c_22[1, 2] /
                Sigma_xxtilde.c_22[2, 2] * (xtilde_i[2] - mu_x.c_i)
              sigsq_xtilde1.xtilde2c <- Sigma_xxtilde.c_22[1, 1] -
                Sigma_xxtilde.c_22[1, 2]^2 / Sigma_xxtilde.c_22[2, 2]
              part2 <- sum(dnorm(x = xtilde_i, log = TRUE,
                                 mean = c(mu_xtilde1.xtilde2c, mu_x.c_i),
                                 sd = sqrt(c(sigsq_xtilde1.xtilde2c,
                                             Sigma_xxtilde.c_22[2, 2]))))

            } else {

              part2 <- dmvnorm(x = xtilde_i, log = TRUE,
                               mean = Mu_xxtilde.c[-1],
                               sigma = Sigma_xxtilde.c_22)

            }

            # Log-likelihood
            ll.vals[ii] <- part1 + part2

          }
          ll.r <- sum(ll.vals)

        } else {

          # Full likelihood

          # Function for integrating out X's
          f_yxtildex.c1 <- function(k_i, y_i, x_i, onec_i, q_i, mu_x.c_i, xtilde_i) {

            x_i <- matrix(x_i, nrow = 1)
            dens <- apply(x_i, 2, function(z) {

              # Transformation
              s_i <- z / (1 - z^2)

              # P(Y|X,C)
              p_y.xc <-
                (1 + exp(-as.numeric(onec_i %*% f.betas[-2, , drop = FALSE]) -
                           s_i * f.beta_x - q_i))^(-1)

              # f(Y,X,Xtilde|C) = f(Y|X,C) f(Xtilde1|Xtilde2,X) f(Xtilde2|X) f(X|C)
              dbinom(x = y_i, size = 1, prob = p_y.xc) *
                prod(dnorm(x = xtilde_i, mean = s_i, sd = sqrt(f.sigsq_m))) *
                dnorm(x = s_i, mean = mu_x.c_i, sd = sqrt(f.sigsq_x.c))

            })

            # Back-transformation
            out <- matrix(dens * (1 + x_i^2) / (1 - x_i^2)^2, ncol = ncol(x_i))
            return(out)

          }

          # Get integration tolerance
          if (estimating.hessian) {
            int.tol <- integrate_tol_hessian
          } else if (identical(f.theta, extra.args$start, ignore.environment = TRUE)) {
            int.tol <- integrate_tol_start
          } else {
            int.tol <- integrate_tol
          }

          int.vals <- c()
          for (ii in 1: length(xtilde.r)) {

            # Get values for ith participant
            k_i <- k.r[ii]
            y_i <- y.r[ii]
            onec_i <- onec.r[ii, ]
            q_i <- q.r[ii]
            mu_x.c_i <- mu_x.c[ii]
            xtilde_i <- xtilde.r[[ii]]

            # Try integrating out X with default settings
            int.ii <-
              adaptIntegrate(f = f_yxtildex.c1, tol = int.tol,
                             lowerLimit = -1, upperLimit = 1,
                             vectorInterface = TRUE,
                             k_i = k_i, y_i = y_i, onec_i = onec_i, q_i = qg_i,
                             mu_x.c_i = mu_x.c_i, xtilde_i = xtilde_i)

            # If integral 0 and f.sigsq_m small, look at region around Xtilde
            if (int.ii$integral == 0 & inside(f.sigsq_m, c(0, 0.1), FALSE)) {

              center.s <- mean(xtilde_i)
              center.x <- (sqrt(4 * center.s^2 + 1) - 1) / (2 * center.s)
              incr <- 1
              for (jj in 1: 6) {
                incr <- incr / 10
                lowupp.x <- c(max(center.x - incr, -1), min(center.x + incr, 1))
                int.ii <-
                  adaptIntegrate(f = f_yxtildex.c1, tol = int.tol,
                                 lowerLimit = lowupp.x[1], upperLimit = lowupp.x[2],
                                 vectorInterface = TRUE,
                                 k_i = k_i, y_i = y_i, onec_i = onec_i, q_i = q_i,
                                 mu_x.c_i = mu_x.c_i, xtilde_i = xtilde_i)
                if (int.ii$integral > 0) {
                  break
                }
              }

            }

            # If integral 0 and f.sigsq_x.c small, look at region around E(X|C)
            if (int.ii$integral == 0 & f.sigsq_x.c < 0.1) {

              center.s <- mu_x.c_i
              center.x <- (sqrt(4 * center.s^2 + 1) - 1) / (2 * center.s)
              incr <- 1
              for (jj in 1: 6) {
                incr <- incr / 10
                lowupp.x <- c(max(center.x - incr, -1), min(center.x + incr, 1))
                int.ii <-
                  adaptIntegrate(f = f_yxtildex.c1, tol = int.tol,
                                 lowerLimit = lowupp.x[1], upperLimit = lowupp.x[2],
                                 vectorInterface = TRUE,
                                 k_i = k_i, y_i = y_i, onec_i = onec_i, q_i = q_i,
                                 mu_x.c_i = mu_x.c_i, xtilde_i = xtilde_i)
                if (int.ii$integral > 0) {
                  break
                }
              }

            }

            int.vals[ii] <- int.ii$integral

            # If integral 0, set skip.rest to TRUE to skip further LL calculations
            if (int.ii$integral == 0) {
              print(paste("Integral is 0 for ii = ", ii, sep = ""))
              print(f.theta)
              print(int.ii)
              skip.rest <- TRUE
              break
            }

          }
          ll.r <- sum(log(int.vals))

        }

      } else {
        ll.r <- 0
      }

      if (! skip.rest & n > 0) {

        # E(X|C)
        mu_x.c <- onec %*% f.alphas

        if (approx_integral) {

          # Probit approximation for logistic-normal integral

          # E(X,Xtilde|C) and V(X,Xtilde|C)
          Mu_xxtilde.c_1 <- mu_x.c
          Mu_xxtilde.c_2 <- mu_x.c
          Sigma_xxtilde.c_11 <- f.sigsq_x.c
          Sigma_xxtilde.c_12 <- f.sigsq_x.c
          Sigma_xxtilde.c_22 <- f.sigsq_x.c + f.sigsq_m

          # E(X|Xtilde,C) and V(X|Xtilde,C)
          mu_x.xtildec <- Mu_xxtilde.c_1 + Sigma_xxtilde.c_12 /
            Sigma_xxtilde.c_22 * (xtilde - Mu_xxtilde.c_2)
          sigsq_x.xtildec <- Sigma_xxtilde.c_11 - Sigma_xxtilde.c_12^2 /
            Sigma_xxtilde.c_22

          # Approximation of \int_x f(Y|X,C) f(X|Xtilde,C) dx
          if (some.cs) {
            t <- as.numeric(f.beta_0 + f.beta_x * mu_x.xtildec +
                              c %*% f.beta_c + q) /
              sqrt(1 + sigsq_x.xtildec * f.beta_x^2 / 1.7^2)
          } else {
            t <- (f.beta_0 + f.beta_x * mu_x.xtildec + q) /
              sqrt(1 + sigsq_x.xtildec * f.beta_x^2 / 1.7^2)
          }
          p <- exp(t) / (1 + exp(t))
          part1 <- dbinom(x = y, size = 1, prob = p, log = TRUE)

          # log[f(Xtilde|C)]
          part2 <- dnorm(x = xtilde, log = TRUE,
                         mean = Mu_xxtilde.c_2, sd = sqrt(Sigma_xxtilde.c_22))

          # Log-likelihood
          ll.vals <- part1 + part2
          ll.s <- sum(ll.vals)

        } else {

          # Full likelihood

          # Function for integrating out X
          f_yxtildex.c2 <- function(y_i, x_i, onec_i, q_i, mu_x.c_i, xtilde_i) {

            # Transformation
            s_i <- x_i / (1 - x_i^2)

            # P(Y|X,C)
            p_y.xc <-
              (1 + exp(as.numeric(-onec_i %*% f.betas[-2, , drop = FALSE]) -
                         s_i * f.beta_x - q_i))^(-1)

            # E(Xtilde|X)
            mu_xtilde.x <- s_i

            # f(Y,X,Xtilde|C) = f(Y|X,C) f(Xtilde|X) f(X|C)
            f_yx.xtildec <-
              dbinom(x = y_i, size = 1, prob = p_y.xc) *
              dnorm(x = xtilde_i, mean = mu_xtilde.x, sd = sqrt(f.sigsq_m)) *
              dnorm(x = s_i, mean = mu_x.c_i, sd = sqrt(f.sigsq_x.c))

            # Back-transformation
            out <- f_yx.xtildec * (1 + x_i^2) / (1 - x_i^2)^2
            return(out)

          }

          # Get integration tolerance
          if (estimating.hessian) {
            int.tol <- integrate_tol_hessian
          } else if (identical(f.theta, extra.args$start, ignore.environment = TRUE)) {
            int.tol <- integrate_tol_start
          } else {
            int.tol <- integrate_tol
          }

          int.vals <- c()
          for (ii in 1: length(xtilde.i)) {

            # Get values for ith participant
            y_i <- y[ii]
            onec_i <- onec[ii, ]
            q_i <- q[ii]
            mu_x.c_i <- mu_x.c[ii]
            xtilde_i <- xtilde.i[ii]

            # Try integrating out X_i with default settings
            int.ii <-
              adaptIntegrate(f = f_yxtildex.c2, tol = int.tol,
                             lowerLimit = -1, upperLimit = 1,
                             vectorInterface = TRUE,
                             y_i = y_i, onec_i = onec_i, q_i = q_i,
                             mu_x.c_i = mu_x.c_i, xtilde_i = xtilde_i)

            # If integral 0 and f.sigsq_m small, look at region around Xtilde
            if (int.ii$integral == 0 & inside(f.sigsq_m, c(0, 0.1), FALSE)) {

              center.s <- xtilde_i
              center.x <- (sqrt(4 * center.s^2 + 1) - 1) / (2 * center.s)
              incr <- 1
              for (jj in 1: 6) {
                incr <- incr / 10
                lowupp.x <- c(max(center.x - incr, -1), min(center.x + incr, 1))
                int.ii <-
                  adaptIntegrate(f = f_yxtildex.c2, tol = int.tol,
                                 lowerLimit = lowupp.x[1], upperLimit = lowupp.x[2],
                                 vectorInterface = TRUE,
                                 y_i = y_i, onec_i = onec_i, q_i = q_i,
                                 mu_x.c_i = mu_x.c_i, xtilde_i = xtilde_i)
                if (int.ii$integral > 0) {
                  break
                }
              }

            }

            # If integral 0 and f.sigsq_x.c small, look at region around E(X|C)
            if (int.ii$integral == 0 & inside(f.sigsq_x.c, c(0, 0.1), FALSE)) {

              center.s <- mu_x.c_i
              center.x <- (sqrt(4 * center.s^2 + 1) - 1) / (2 * center.s)
              incr <- 1
              for (jj in 1: 6) {
                incr <- incr / 10
                lowupp.x <- c(max(center.x - incr, -1), min(center.x + incr, 1))
                int.ii <-
                  adaptIntegrate(f = f_yxtildex.c2, tol = int.tol,
                                 lowerLimit = lowupp.x[1], upperLimit = lowupp.x[2],
                                 vectorInterface = TRUE,
                                 y_i = y_i, onec_i = onec_i, q_i = q_i,
                                 mu_x.c_i = mu_x.c_i, xtilde_i = xtilde_i)
                if (int.ii$integral > 0) {
                  break
                }
              }

            }

            int.vals[ii] <- int.ii$integral

            # If integral 0, set skip.rest to TRUE to skip further LL calculations
            if (int.ii$integral == 0) {
              print(paste("Integral is 0 for ii = ", ii, sep = ""))
              print(f.theta)
              print(int.ii)
              skip.rest <- TRUE
              break
            }

          }
          ll.s <- sum(log(int.vals))

        }

      }

      ll <- ll.r + ll.s

    } else {

      # Likelihood:
      # L = f(Y|X,C) f(X|C)

      # P(Y|X,C)
      eta <- onexc %*% f.betas
      p_y.xc <- (1 + exp(-eta))^(-1)

      # E(X|C)
      mu_x.c <- onec %*% f.alphas

      # Log-likelihood
      ll <- sum(dbinom(x = y, log = TRUE,
                       size = 1, prob = p_y.xc) +
                  dnorm(x = xtilde, log = TRUE, mean = mu_x.c, sd = sqrt(f.sigsq_x.c)))

    }

    # Return negative log-likelihood
    return(-ll)

  }

  # Create list of extra arguments, and assign default starting values and
  # lower values if not specified by user
  extra.args <- list(...)
  if (is.null(extra.args$start)) {
    if (! merror) {
      extra.args$start <- c(rep(0.01, n.betas + n.alphas), 1)
    } else {
      extra.args$start <- c(rep(0.01, n.betas + n.alphas), 1, 1)
    }
  }
  if (is.null(extra.args$lower)) {
    if (! merror) {
      extra.args$lower <- c(rep(-Inf, n.betas + n.alphas), 1e-3)
    } else {
      extra.args$lower <- c(rep(-Inf, n.betas + n.alphas), rep(1e-3, 2))
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

  # Create list to return
  theta.hat <- ml.max$par
  names(theta.hat) <- theta.labels
  ret.list <- list(theta.hat = theta.hat)

  # If requested, add variance-covariance matrix to ret.list
  if (estimate_var) {
    hessian.mat <- pracma::hessian(f = ll.f, estimating.hessian = TRUE,
                                   x0 = theta.hat)
    theta.variance <- try(solve(hessian.mat), silent = TRUE)
    if (class(theta.variance) == "try-error") {
      print(hessian.mat)
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


# # Prospective sampling
# n <- 10000
# n.reps <- n / 10
# betas <- c(-2, 0.5)
# alphas <- 0
# sigsq_x.c <- 1
# sigsq_m <- 0.25
#
# prev <- samp_y1y0 <- NULL
# merror <- TRUE
#
# x <- rnorm(n, mean = alphas, sd = sqrt(sigsq_x.c))
# y <- rbinom(n = n, size = 1, prob = (1 + exp(-betas[1] - betas[2] * x))^(-1))
#
# xtilde <- list()
# locs.2 <- sample(1: n, n.reps, replace = FALSE)
# for (ii in 1: length(x)) {
#   xtilde[[ii]] <- x[ii] + rnorm(n = ifelse(ii %in% locs.rep, 2, 1))
# }
#
# abc <- logreg_xerrors(y = y, xtilde = xtilde, control = list(trace = 1))
#
# # Try to get X very different in cases and controls
# n <- 1000000
# betas <- c(-10, 5)
# alphas <- 0
# sigsq_x.c <- 1
# x <- rnorm(n, mean = alphas, sd = sqrt(sigsq_x.c))
# y <- rbinom(n = n, size = 1, prob = (1 + exp(-betas[1] - betas[2] * x))^(-1))
# mean(y)
# par(mfrow = c(3, 1))
# histo(x[y == 1], xlim = c(-4, 4), breaks = 25, dis = "norm")
# histo(x[y == 0], xlim = c(-4, 4), breaks = 25, dis = "norm")
# histo(c(x[which(y == 1)[1: 1000]], x[which(y == 0)[1: 1000]]), breaks = 25, dis = "norm")
#
# # Case-control sampling
# n.pop <- 5000000
# n.samp <- 50000
# n.reps <- n.samp / 20
# betas <- c(-10, 5)
# alphas <- 0
# sigsq_x.c <- 1
# sigsq_m <- 0.1
#
# x <- rnorm(n = n.pop, mean = alphas, sd = sqrt(sigsq_x.c))
# y <- rbinom(n = n.pop, size = 1, prob = (1 + exp(-betas[1] - betas[2] * x))^(-1))
#
# loc.cases <- which(y == 1)
# loc.controls <- which(y == 0)
# loc.sampled <- c(sample(loc.cases, n.samp, replace = FALSE), sample(loc.controls, n.samp, replace = FALSE))
# x <- x[loc.sampled]
# y <- y[loc.sampled]
#
# xtilde <- list()
# locs.2 <- sample(1: (n.samp * 2), n.reps * 2, replace = FALSE)
# for (ii in 1: length(x)) {
#   xtilde[[ii]] <- x[ii] + rnorm(n = ifelse(ii %in% locs.2, 2, 1))
# }
#
# summary(glm(y ~ x, family = "binomial"))
# summary(glm(y ~ x, offset = rep(log(0.97/0.03), length(x)), family = "binomial"))
# summary(glm(y ~ sapply(xtilde, function(x) x[1]), offset = rep(log(0.97/0.03), length(x)), family = "binomial"))
# abc <- logreg_xerrors(y = y, xtilde = xtilde, prev = 0.03, control = list(trace = 1))
#
# x1 <- sapply(xtilde, function(x) x[1])
# x2 <- sapply(xtilde, function(x) x[2])
# plot(x1, x2)
#
# # Back to prospective sampling for this new scenario
# n <- 10000
# n.reps <- n / 10
# betas <- c(-10, 5)
# alphas <- 0
# sigsq_x.c <- 1
# sigsq_m <- 0.1
#
# x <- rnorm(n, mean = alphas, sd = sqrt(sigsq_x.c))
# y <- rbinom(n = n, size = 1, prob = (1 + exp(-betas[1] - betas[2] * x))^(-1))
#
# xtilde <- list()
# locs.2 <- sample(1: n, n.reps, replace = FALSE)
# for (ii in 1: length(x)) {
#   xtilde[[ii]] <- x[ii] + rnorm(n = ifelse(ii %in% locs.rep, 2, 1))
# }
#
# abc <- logreg_xerrors(y = y, xtilde = xtilde, control = list(trace = 1))
#
# # Loop through sampling rates and record estimates for single large-n trial
# betas <- c(-10, 5)
# alphas <- 0
# sigsq_x.c <- 1
# sigsq_m <- 0.1
#
# n.pop <- 5000000
# n.samp <- 100000
# n.samp <- 100000
# n.reps <- n.samp / 50
# n.reps <- 500
# vals <- c(seq(0.01, 0.1, 0.01), seq(0.2, 0.9, 0.1))
#
# estimates <- matrix(NA, nrow = length(vals), ncol = 5)
# for (ii in 1: length(vals)) {
#
#   p_case <- vals[ii]
#
#   x <- rnorm(n = n.pop, mean = alphas, sd = sqrt(sigsq_x.c))
#   y <- rbinom(n = n.pop, size = 1, prob = (1 + exp(-betas[1] - betas[2] * x))^(-1))
#
#   loc.cases <- which(y == 1)
#   loc.controls <- which(y == 0)
#   loc.sampled <- c(sample(loc.cases, n.samp * p_case, replace = FALSE),
#                    sample(loc.controls, n.samp * (1 - p_case), replace = FALSE))
#   x <- x[loc.sampled]
#   y <- y[loc.sampled]
#
#   xtilde <- list()
#   locs.2 <- sample(1: n.samp, n.reps, replace = FALSE)
#   for (jj in 1: length(x)) {
#     xtilde[[jj]] <- x[jj] + rnorm(n = ifelse(jj %in% locs.2, 2, 1), sd = sqrt(sigsq_m))
#   }
#
#   #summary(glm(y ~ x, family = "binomial"))
#   #summary(glm(y ~ x, offset = rep(log(0.97/0.03), length(x)), family = "binomial"))
#   #summary(glm(y ~ sapply(xtilde, function(x) x[1]), offset = rep(log(0.97/0.03), length(x)), family = "binomial"))
#   abc <- logreg_xerrors(y = y, xtilde = xtilde, prev = 0.03, control = list(trace = 1))
#   estimates[ii, ] <- abc$theta.hat
#
# }
