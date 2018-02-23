#' Poolwise Logistic Regression with Gamma-Distributed Exposure Subject to
#' Errors
#'
#' Similar to \code{\link{p.logreg.xerrors}}, but for skewed exposure that takes
#' on only positive values. Assumes a constant-scale Gamma model for exposure
#' given covariates, and multiplicative lognormal processing errors and
#' measurement errors acting on the poolwise mean exposure. Borrows ideas from
#' [1-3]. Manuscript fully describing the approach is in preparation.
#'
#'
#' @param g Numeric vector of pool sizes (number of individuals in each pool).
#' Do not have to specify if \code{c} is specified (function can figure it out).
#'
#' @param y Numeric vector with poolwise \code{Y} values, coded 0 if all members
#' are controls and 1 if all members are cases.
#'
#' @param xtilde Numeric vector (or list of numeric vectors, if some pools have
#' replicate measurements) with poolwise \code{Xtilde} values.
#'
#' @param c List where each element is a matrix containing the \code{\bold{C}}
#' values for individual members of a particular pool. Each row of the matrix
#' should contain the \code{\bold{C}} values for one individual. Can leave as
#' \code{NULL} if there are no \code{\bold{C}} variables.
#'
#' @param error.type Character string specifying the errors that \code{X} is
#' subject to.  Possible values are \code{"neither"} for neither processing
#' error nor measurement error, \code{"processing"} for processing error only,
#' \code{"measurement"} for measurement error only, and \code{"both"} for both.
#'
#' @param p_y1 Optional numeric value specifying disease prevalence, allowing
#' for valid estimation of the intercept with case-control sampling.  If it's
#' easier, you can specify \code{p_sample_y1y0} instead.  Only used if
#' \code{offset.formula = 1}.
#'
#' @param p_sample.y1y0 Optional numeric vector if length 2 specifying sampling
#' probabilities for cases and controls, allowing for valid estimation of the
#' intercept with case-control sampling.  If it's easier, you can specify
#' \code{p_y1} instead.
#'
#' @param integrate.tol Numeric value specifying the \code{tol} input to
#' \code{\link{adaptIntegrate}}.  Only used if \code{approx.integral = FALSE}.
#'
#' @param integrate.tol.start Same as \code{integrate.tol}, but applies only to
#' the very first iteration of ML maximization.  The first iteration tends to
#' take much longer than subsequent ones, so less precise integration at the
#' start can speed things up.
#'
#' @param integrate.tol.hessian Same as \code{integrate.tol}, but for use when
#' estimating the Hessian matrix only.  Sometimes more precise integration
#' (i.e. smaller tolerance) than used for maximizing the likelihood helps
#' prevent cases where the inverse Hessian is not positive definite.
#'
#' @param estimate.var If \code{TRUE}, function returns variance-covariance
#' matrix for parameter estimates, calculated as the inverse of the estimated
#' Hessian matrix at the MLE's.
#'
#' #' @return A list containing the following:
#' \enumerate{
#'   \item Numeric vector of parameter estimates.
#'   \item Variance-covariance matrix (if \code{estimate.var = TRUE}).
#'   \item The returned \code{\link[stats]{nlminb}} object from maximizing the
#' log-likelihood function.
#'   \item Akaike information criterion (AIC).
#' }
#'
#'
#' @export
p_logreg_xerrors_gamma <- function(g = NULL, y, xtilde, c = NULL,
                                   errors = "both",
                                   diff.pe = FALSE, diff.me = FALSE,
                                   p_y1 = NULL, p_sample.y1y0 = NULL,
                                   integrate.tol = 1e-8,
                                   integrate.tol.start = integrate.tol,
                                   integrate.tol.hessian = integrate.tol,
                                   estimate.var = FALSE, ...) {

  # Check that inputs are valid
  if (! error.type %in% c("neither", "processing", "measurement", "both")) {
    stop("The input 'error.type' should be set to 'neither', 'processing',
         'measurement', or 'both'.")
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
  if (! (is.numeric(integrate.tol) & inside(integrate.tol, c(1e-32, Inf)))) {
    stop("The input 'integrate.tol' must be a numeric value greater than 1e-32.")
  }
  if (! (is.numeric(integrate.tol.start) & inside(integrate.tol.start, c(1e-32, Inf)))) {
    stop("The input 'integrate.tol.start' must be a numeric value greater than 1e-32.")
  }
  if (! (is.numeric(integrate.tol.hessian) & inside(integrate.tol.hessian, c(1e-32, Inf)))) {
    stop("The input 'integrate.tol.hessian' must be a numeric value greater than 1e-32.")
  }
  if (! is.logical(estimate.var)) {
    stop("The input 'estimate.var' should be TRUE or FALSE.")
  }

  # Get name of xtilde input
  x.varname <- deparse(substitute(xtilde))

  # Get information about covariates C
  if (is.null(c)) {
    c.varnames <- NULL
    n.cvars <- 0
    some.cs <- FALSE
  } else {
    if (! is.matrix(c)) {
      c <- as.matrix(c)
    }
    n.cvars <- ncol(c[[1]])
    some.cs <- TRUE
    c.varnames <- colnames(c[[1]])
    if (is.null(c.varnames)) {
      if (n.cvars == 1) {
        c.varnames <- deparse(substitute(c))
      } else {
        c.varnames <- paste("c", 1: n.cvars, sep = "")
      }
    }
  }

  # Get number of betas and lambdas
  n.betas <- 2 + n.cvars
  n.lambdas <- 1 + n.cvars

  # Figure out pool sizes if not specified
  if (is.null(g)) {
    g <- sapply(c, nrow)
  }

  # Create vector indicating which observations are pools
  ipool <- ifelse(g > 1, 1, 0)

  # Construct X|C design matrix list
  if (some.cs) {
    design <- lapply(c, function(x) cbind(1, x))
  } else {
    design <- NULL
  }

  # Calculate poolwise C's
  if (some.cs) {
    ci <- matrix(sapply(c, colSums), byrow = TRUE, ncol = n.cvars)
  } else {
    ci <- NULL
  }

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

  # Separate out pools with precisely measured X
  if (error.type == "neither") {
    which.p <- 1: n.pools
  } else if (error.type == "processing") {
    which.p <- which(ipool == 0)
  } else {
    which.p <- NULL
  }
  n.p <- length(which.p)
  some.p <- n.p > 0
  if (some.p) {
    g.p <- g[which.p]
    y.p <- y[which.p]
    x.p <- unlist(xtilde[which.p])
    ci.p <- ci[which.p, , drop = FALSE]
    gxc.p <- cbind(g.p, x.p, ci.p)
    qg.p <- qg[which.p]
    design.p <- design[which.p]
  }

  # Separate out pools with replicate Xtilde measurements
  class.xtilde <- class(xtilde)
  if (class.xtilde == "list") {
    k <- sapply(xtilde, length)
    which.r <- which(k > 1)
    n.r <- length(which.r)
    some.r <- n.r > 0
    if (some.r) {
      k.r <- k[which.r]
      g.r <- g[which.r]
      ipool.r <- ipool[which.r]
      y.r <- y[which.r]
      xtilde.r <- xtilde[which.r]
      design.r <- design[which.r]
      ci.r <- ci[which.r, , drop = FALSE]
      qg.r <- qg[which.r]
    }
  } else {
    k <- rep(1, n.pools)
    some.r <- FALSE
  }

  # Separate out pools with single imprecisely measured X
  if (error.type == "neither") {
    which.i <- NULL
  } else if (error.type == "processing") {
    which.i <- which(ipool == 1 & k == 1)
  } else if (error.type %in% c("measurement", "both")) {
    which.i <- which(k == 1)
  }
  n.i <- length(which.i)
  some.i <- n.i > 0
  if (some.i) {
    g.i <- g[which.i]
    ipool.i <- ipool[which.i]
    y.i <- y[which.i]
    ci.i <- ci[which.i, , drop = FALSE]
    qg.i <- qg[which.i]
    design.i <- design[which.i]
    xtilde.i <- unlist(xtilde[which.i])
  }

  # Get indices for parameters being estimated, and create labels as well
  loc.betas <- 1: n.betas
  beta.labels <- paste("beta", c("0", x.varname, c.varnames), sep = "_")

  loc.lambdas <- (n.betas + 1): (n.betas + n.lambdas)
  lambda.labels <- paste("lambda", c("0", c.varnames), sep = "_")

  loc.b <- n.betas + n.lambdas + 1

  if (error.type == "neither") {
    theta.labels <- c(beta.labels, lambda.labels, "b")
  } else if (error.type == "processing") {
    if (diff.pe) {
      loc.sigsq_p1 <- loc.b + 1
      loc.sigsq_p0 <- loc.b + 2
      theta.labels <- c(beta.labels, lambda.labels,
                        "b", "sigsq_p1", "sigsq_p0")
    } else {
      loc.sigsq_p1 <- loc.sigsq_p0 <- loc.b + 1
      theta.labels <- c(beta.labels, lambda.labels, "b", "sigsq_p")
    }
  } else if (error.type == "measurement") {
    if (diff.me) {
      loc.sigsq_m1 <- loc.b + 1
      loc.sigsq_m0 <- loc.b + 2
      theta.labels <- c(beta.labels, lambda.labels,
                        "b", "sigsq_m1", "sigsq_m0")
    } else {
      loc.sigsq_m1 <- loc.sigsq_m0 <- loc.b + 1
      theta.labels <- c(beta.labels, lambda.labels, "b", "sigsq_m")
    }
  } else if (error.type == "both") {
    if (diff.pe & diff.me) {
      loc.sigsq_p1 <- loc.b + 1
      loc.sigsq_p0 <- loc.b + 2
      loc.sigsq_m1 <- loc.b + 3
      loc.sigsq_m0 <- loc.b + 4
      theta.labels <- c(beta.labels, lambda.labels,
                        "b", "sigsq_p1", "sigsq_p0", "sigsq_m1",
                        "sigsq_m0")
    } else if (diff.pe & ! diff.me) {
      loc.sigsq_p1 <- loc.b + 1
      loc.sigsq_p0 <- loc.b + 2
      loc.sigsq_m1 <- loc.sigsq_m0 <- loc.b + 3
      theta.labels <- c(beta.labels, lambda.labels,
                        "b", "sigsq_p1", "sigsq_p0", "sigsq_m")
    } else if (! diff.pe & diff.me) {
      loc.sigsq_p1 <- loc.sigsq_p0 <- loc.b + 1
      loc.sigsq_m1 <- loc.b + 2
      loc.sigsq_m0 <- loc.b + 3
      theta.labels <- c(beta.labels, lambda.labels,
                        "b", "sigsq_p", "sigsq_m1", "sigsq_m0")
    } else if (! diff.pe & ! diff.me) {
      loc.sigsq_p1 <- loc.sigsq_p0 <- loc.b + 1
      loc.sigsq_m1 <- loc.sigsq_m0 <- loc.b + 2
      theta.labels <- c(beta.labels, lambda.labels,
                        "b", "sigsq_p", "sigsq_m")
    }
  }

  # Log-likelihood function
  ll.f <- function(f.theta, estimating.hessian = FALSE) {

    # Extract parameters
    f.betas <- matrix(f.theta[loc.betas], ncol = 1)
    f.beta_0 <- f.betas[1]
    f.beta_x <- f.betas[2]
    f.beta_c <- matrix(f.betas[-c(1: 2)], ncol = 1)

    f.lambdas <- matrix(f.theta[loc.lambdas], ncol = 1)
    f.lambda_0 <- f.lambdas[1]
    f.lambda_c <- matrix(f.lambdas[-1], ncol = 1)

    f.b <- f.theta[loc.b]

    if (error.type == "neither") {
      f.sigsq_p1 <- f.sigsq_p0 <- f.sigsq_m1 <- f.sigsq_m0 <- 0
    }
    if (error.type %in% c("processing", "both")) {
      f.sigsq_p1 <- f.theta[loc.sigsq_p1]
      f.sigsq_p0 <- f.theta[loc.sigsq_p0]
    } else {
      f.sigsq_p1 <- f.sigsq_p0 <- 0
    }
    if (error.type %in% c("measurement", "both")) {
      f.sigsq_m1 <- f.theta[loc.sigsq_m1]
      f.sigsq_m0 <- f.theta[loc.sigsq_m0]
    } else {
      f.sigsq_m1 <- f.sigsq_m0 <- 0
    }

    if (some.p) {

      # Likelihood for pools with precisely measured X_i^*:
      # L = f(Y_i|X_i^*,C_i^*) f(X_i^*|C_i1, ..., C_igi)

      # P(Y_i|X_i^*,C_i^*)
      eta <- gxc.p %*% f.betas + qg.p
      p_y.xc <- (1 + exp(-eta))^(-1)

      # a_i's in X_|C_i1, ..., C_igi ~ Gamma(a_i, b)
      if (some.cs) {
        alphas <- sapply(design.p, function(x) sum(exp(x %*% f.lambdas)))
      } else {
        alphas <- g.p * exp(f.lambda_0)
      }

      # Log-likelihood (non-positive X that generate NAs get excluded)
      ll.p <- sum(log(dbinom(x = y.p, size = 1, prob = p_y.xc) *
                        dgamma(x = x.p, shape = alphas, scale = f.b)))

    } else {
      ll.p <- 0
    }

    # Set skip.rest flag to FALSE
    skip.rest <- FALSE

    if (some.r) {

      # Likelihood for pools with replicate Xtilde_i^*'s:
      # L_i = \int_X_i^* f(Y_i|X_i^*,C_i^*) f(Xtilde_i^*|X_i^*)
      #                  f(X_i^*|C_i1, ..., C_igi) dX_i^*

      # Create error vectors
      sigsq_p <- ifelse(y.r > 0, f.sigsq_p1, f.sigsq_p0) * ipool.r
      sigsq_m <- ifelse(y.r > 0, f.sigsq_m1, f.sigsq_m0)

      # a_i's to feed to integral
      if (some.cs) {
        alphas <- sapply(design.r, function(x) sum(exp(x %*% f.lambdas)))
      } else {
        alphas <- g.r * exp(f.lambda_0)
      }

      # Function for integrating out X_i^*'s
      int.f_i1 <- function(k_i, g_i, ipool_i, y_i, x_i, ci_i, qg_i, xtilde_i,
                           a_i, sigsq_p_i, sigsq_m_i) {

        # f(Y_i,Xtilde_i^*,X_i^*|C_i1, ..., C_ig)
        x_i <- matrix(x_i, nrow = 1)
        f_yxtildex.c <- apply(x_i, 2, function(z) {

          # Transformation
          s_i <- z / (1 - z)

          # P(Y|X,C)
          if (some.cs) {
            p_y.xc <- (1 + exp(-g_i * f.beta_0 - f.beta_x * s_i -
                                 as.vector(t(f.beta_c) %*% ci_i) - qg_i))^(-1)
          } else {
            p_y.xc <- (1 + exp(-g_i * f.beta_0 - f.beta_x * s_i - qg_i))^(-1)
          }

          # E(log(e_i)) and V(log(e_i))
          Mu_loge <- rep(-1/2 * (sigsq_p_i * ipool_i + sigsq_m_i), k_i)
          Sigma_loge <- matrix(sigsq_p_i * ipool_i, ncol = k_i, nrow = k_i) +
            diag(x = sigsq_m_i, ncol = k_i, nrow = k_i)

          # Density
          dbinom(x = y_i, size = 1, prob = p_y.xc) *
            1 / prod(xtilde_i) * dmvnorm(x = log(1 / s_i * xtilde_i),
                                         mean = Mu_loge, sigma = Sigma_loge) *
            dgamma(x = s_i, shape = a_i, scale = f.b)

        })

        # Back-transformation
        out <- matrix(f_yxtildex.c / (1 - x_i)^2, ncol = ncol(x_i))

      }

      # Get integration tolerance
      if (estimating.hessian) {
        int.tol <- integrate.tol.hessian
      } else if (all(f.theta == extra.args$start)) {
        int.tol <- integrate.tol.start
      } else {
        int.tol <- integrate.tol
      }

      int.vals <- c()
      for (ii in 1: length(xtilde.r)) {

        # Get values for ith participant
        g_i <- g.r[ii]
        ipool_i <- ipool.r[ii]
        k_i <- k.r[ii]
        y_i <- y.r[ii]
        ci_i <- ci.r[ii, ]
        qg_i <- qg.r[ii]
        xtilde_i <- xtilde.r[[ii]]
        a_i <- alphas[ii]
        sigsq_p_i <- sigsq_p[ii]
        sigsq_m_i <- sigsq_m[ii]

        # Try integrating out X_i with default settings
        int.ii <-
          adaptIntegrate(f = int.f_i1, tol = int.tol,
                         lowerLimit = 0, upperLimit = 1,
                         vectorInterface = TRUE,
                         g_i = g_i, ipool_i = ipool_i, k_i = k_i, y_i = y_i,
                         ci_i = ci_i, qg_i = qg_i, xtilde_i = xtilde_i,
                         a_i = a_i, sigsq_p_i = sigsq_p_i,
                         sigsq_m_i = sigsq_m_i)
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

    } else {
      ll.r <- 0
    }

    if (some.i & ! skip.rest) {

      # Likelihood for pools with single imprecisely measured Xtilde_i^*:
      # L_i = \int_X_i^* f(Y_i|X_i^*,C_i^*) f(Xtilde_i^*|X_i^*)
      #                  f(X_i^*|C_i1, ..., C_igi) dX_i^*

      # Create error vectors
      sigsq_p <- ifelse(y.i > 0, f.sigsq_p1, f.sigsq_p0) * ipool.i
      sigsq_m <- ifelse(y.i > 0, f.sigsq_m1, f.sigsq_m0)

      # a_i's to feed to integral
      if (some.cs) {
        alphas <- sapply(design.i, function(x) sum(exp(x %*% f.lambdas)))
      } else {
        alphas <- g.i * exp(f.lambda_0)
      }

      # Function for integrating out X_i^*'s
      int.f_i2 <- function(g_i, ipool_i, y_i, x_i, ci_i, qg_i, xtilde_i, a_i,
                           sigsq_p_i, sigsq_m_i) {

        # Transformation
        s_i <- x_i / (1 - x_i)

        # P(Y_i|X_i^*,C_i^*)
        if (some.cs) {
          p_y.xc <- (1 + exp(-g_i * f.beta_0 - f.beta_x * s_i -
                               as.vector(t(f.beta_c) %*% ci_i) - qg_i))^(-1)
        } else {
          p_y.xc <- (1 + exp(-g_i * f.beta_0 - f.beta_x * s_i - qg_i))^(-1)
        }

        # E(log(e_i)) and V(log(e_i))
        mu_loge <- -1/2 * (sigsq_p_i * ipool_i + sigsq_m_i)
        sigsq_loge <- sigsq_p_i * ipool_i + sigsq_m_i

        # f(Y_i,Xtilde_i^*,X_i^*|C_i1, ..., C_igi) = f(Y_i|X_i^*,C_i^*)
        # f(Xtilde_i^*|X_i^*) f(X_i^*|C_i1, ..., C_igi)
        f_yxtildex.c <- dbinom(x = y_i, size = 1, prob = p_y.xc) *
          1 / xtilde_i * dnorm(x = log(xtilde_i / s_i),
                               mean = mu_loge, sd = sqrt(sigsq_loge)) *
          dgamma(x = s_i, shape = a_i, scale = f.b)

        # Back-transformation
        out <- f_yxtildex.c / (1 - x_i)^2
        return(out)

      }

      # Get integration tolerance
      if (estimating.hessian) {
        int.tol <- integrate.tol.hessian
      } else if (all(f.theta == extra.args$start)) {
        int.tol <- integrate.tol.start
      } else {
        int.tol <- integrate.tol
      }

      int.vals <- c()
      for (ii in 1: length(xtilde.i)) {

        # Get values for ith participant
        g_i <- g.i[ii]
        ipool_i <- ipool.i[ii]
        y_i <- y.i[ii]
        ci_i <- ci.i[ii, ]
        qg_i <- qg.i[ii]
        xtilde_i <- xtilde.i[ii]
        a_i <- alphas[ii]
        sigsq_p_i <- sigsq_p[ii]
        sigsq_m_i <- sigsq_m[ii]

        # Try integrating out X_i with default settings
        int.ii <-
          adaptIntegrate(f = int.f_i2, tol = int.tol,
                         lowerLimit = 0, upperLimit = 1,
                         vectorInterface = TRUE,
                         g_i = g_i, ipool_i = ipool_i, y_i = y_i,
                         ci_i = ci_i, qg_i = qg_i, xtilde_i = xtilde_i,
                         a_i = a_i, sigsq_p_i = sigsq_p_i,
                         sigsq_m_i = sigsq_m_i)
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
      ll.i <- sum(log(int.vals))

    } else {
      ll.i <- 0
    }

    # Return negative log-likelihood
    ll <- ll.p + ll.r + ll.i
    return(-ll)

  }

  # Create list of extra arguments, and assign default starting values and lower
  # values if not specified by user
  extra.args <- list(...)
  if (is.null(extra.args$start)) {
    if (error.type == "neither") {
      extra.args$start <- c(rep(0, n.betas + n.lambdas), 1)
    } else if (error.type == "processing") {
      extra.args$start <- c(rep(0, n.betas + n.lambdas),
                            rep(1, loc.sigsq_p0 - loc.b + 1))
    } else if (error.type %in% c("measurement", "both")) {
      extra.args$start <- c(rep(0, n.betas + n.lambdas),
                            rep(1, loc.sigsq_m0 - loc.b + 1))
    }
  }
  if (is.null(extra.args$lower)) {
    if (error.type == "neither") {
      extra.args$lower <- c(rep(-Inf, n.betas + n.lambdas), 1e-3)
    } else if (error.type == "processing") {
      extra.args$lower <- c(rep(-Inf, n.betas + n.lambdas),
                            rep(1e-3, loc.sigsq_p0 - loc.b + 1))
    } else if (error.type %in% c("measurement", "both")) {
      extra.args$lower <- c(rep(-Inf, n.betas + n.lambdas),
                            rep(1e-3, loc.sigsq_m0 - loc.b + 1))
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
  if (estimate.var) {
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

  # return ret.list
  return(ret.list)

}
