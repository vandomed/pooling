#' Poolwise Logistic Regression with Gamma Exposure Subject to Errors
#'
#' Assumes constant-scale Gamma model for exposure given covariates, and
#' multiplicative lognormal processing errors and measurement errors acting on
#' the poolwise mean exposure. Manuscript fully describing the approach is
#' under review.
#'
#'
#' @inheritParams p_logreg_xerrors
#'
#' @param c List where each element is a numeric matrix containing the
#' \strong{\code{C}} values for members of a particular pool (1 row for each
#' member).
#' @param integrate_tol Numeric value specifying the \code{tol} input to
#' \code{\link{hcubature}}.
#'
#'
#' @inherit p_logreg_xerrors return
#'
#'
#' @references
#' Mitchell, E.M, Lyles, R.H., and Schisterman, E.F. (2015) "Positing, fitting,
#' and selecting regression models for pooled biomarker data." \emph{Stat. Med}
#' \strong{34}(17): 2544--2558.
#'
#' Schisterman, E.F., Vexler, A., Mumford, S.L. and Perkins, N.J. (2010) "Hybrid
#' pooled-unpooled design for cost-efficient measurement of biomarkers."
#' \emph{Stat. Med.} \strong{29}(5): 597--613.
#'
#' Weinberg, C.R. and Umbach, D.M. (1999) "Using pooled exposure assessment to
#' improve efficiency in case-control studies." \emph{Biometrics} \strong{55}:
#' 718--726.
#'
#' Weinberg, C.R. and Umbach, D.M. (2014) "Correction to 'Using pooled exposure
#' assessment to improve efficiency in case-control studies' by Clarice R.
#' Weinberg and David M. Umbach; 55, 718--726, September 1999."
#' \emph{Biometrics} \strong{70}: 1061.
#'
#' Whitcomb, B.W., Perkins, N.J., Zhang, Z., Ye, A., and Lyles, R. H. (2012)
#' "Assessment of skewed exposure in case-control studies with pooling."
#' \emph{Stat. Med.} \strong{31}: 2461--2472.
#'
#'
#' @examples
#' # Load dataset with (g, Y, Xtilde, C) values for 248 pools and list of C
#' # values for members of each pool. Xtilde values are affected by processing
#' # error.
#' data(pdat2)
#' dat <- pdat2$dat
#' c.list <- pdat2$c.list
#'
#' # Estimate log-OR for X and Y adjusted for C, ignoring processing error
#' fit1 <- p_logreg_xerrors2(
#'   g = dat$g,
#'   y = dat$y,
#'   xtilde = dat$xtilde,
#'   c = c.list,
#'   errors = "neither"
#' )
#' fit1$theta.hat
#'
#' # Repeat, but accounting for processing error.
#' \dontrun{
#' fit2 <- p_logreg_xerrors2(
#'   g = dat$g,
#'   y = dat$y,
#'   xtilde = dat$xtilde,
#'   c = c.list,
#'   errors = "processing",
#'   control = list(trace = 1)
#' )
#' fit2$theta.hat
#' }
#'
#'
#' @export
p_logreg_xerrors2 <- function(g = NULL, y, xtilde, c = NULL,
                              errors = "processing",
                              nondiff_pe = TRUE, nondiff_me = TRUE,
                              constant_pe = TRUE,
                              prev = NULL, samp_y1y0 = NULL,
                              integrate_tol = 1e-8,
                              integrate_tol_start = integrate_tol,
                              integrate_tol_hessian = integrate_tol,
                              estimate_var = TRUE,
                              fix_posdef = FALSE,
                              ...) {

  # Check that inputs are valid
  if (! errors %in% c("neither", "processing", "measurement", "both")) {
    stop("The input 'errors' should be set to 'neither', 'processing',
         'measurement', or 'both'.")
  }
  if (! is.logical(nondiff_pe)) {
    stop("The input 'nondiff_pe' should be TRUE if you want to assume non-differential processing error and FALSE otherwise.")
  }
  if (! is.logical(nondiff_me)) {
    stop("The input 'nondiff_me' should be TRUE if you want to assume non-differential measurement error and FALSE otherwise.")
  }
  if (! is.logical(constant_pe)) {
    stop("The input 'constant_pe' should be TRUE if you want to assume that processing error variance is constant with pool size and FALSE otherwise.")
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
  if (! (is.numeric(integrate_tol) & inside(integrate_tol, c(1e-32, Inf)))) {
    stop("The input 'integrate_tol' must be a numeric value greater than 1e-32.")
  }
  if (! (is.numeric(integrate_tol_start) & inside(integrate_tol_start, c(1e-32, Inf)))) {
    stop("The input 'integrate_tol_start' must be a numeric value greater than 1e-32.")
  }
  if (! (is.numeric(integrate_tol_hessian) & inside(integrate_tol_hessian, c(1e-32, Inf)))) {
    stop("The input 'integrate_tol_hessian' must be a numeric value greater than 1e-32.")
  }
  if (! is.logical(estimate_var)) {
    stop("The input 'estimate_var' should be TRUE or FALSE.")
  }

  # Get name of xtilde input
  x.varname <- deparse(substitute(xtilde))
  if (length(grep("$", x.varname, fixed = TRUE)) > 0) {
    x.varname <- substr(x.varname,
                        start = which(unlist(strsplit(x.varname, "")) == "$") + 1,
                        stop = nchar(x.varname))
  }

  # Get information about covariates C
  if (is.null(c)) {
    c.varnames <- NULL
    n.cvars <- 0
    some.cs <- FALSE
  } else {
    c.varname <- deparse(substitute(c))
    n.cvars <- ncol(c[[1]])
    some.cs <- TRUE
    c.varnames <- colnames(c[[1]])
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

  # Get number of betas and alphas
  n.betas <- 2 + n.cvars
  n.alphas <- 1 + n.cvars

  # Figure out pool sizes if not specified
  if (is.null(g)) {
    g <- sapply(c, nrow)
  }

  # Create vector indicating which observations are pools
  Ig <- ifelse(g > 1, 1, 0)

  # Construct X|C design matrix list
  if (some.cs) {
    onec <- lapply(c, function(x) cbind(1, x))
  } else {
    onec <- NULL
  }

  # Calculate poolwise C's
  if (some.cs) {
    cstar <- matrix(sapply(c, colSums), byrow = TRUE, ncol = n.cvars)
  } else {
    cstar <- NULL
  }

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

  # Separate out pools with precisely measured X
  if (errors == "neither") {
    which.p <- 1: n
  } else if (errors == "processing") {
    which.p <- which(Ig == 0)
  } else {
    which.p <- NULL
  }
  n.p <- length(which.p)
  some.p <- n.p > 0
  if (some.p) {
    g.p <- g[which.p]
    y.p <- y[which.p]
    x.p <- unlist(xtilde[which.p])
    cstar.p <- cstar[which.p, , drop = FALSE]
    gxc.p <- cbind(g.p, x.p, cstar.p)
    qg.p <- qg[which.p]
    onec.p <- onec[which.p]
  }

  # Separate out pools with replicates
  class.xtilde <- class(xtilde)
  if (class.xtilde == "list") {
    k <- sapply(xtilde, length)
    which.r <- which(k > 1)
    n.r <- length(which.r)
    some.r <- n.r > 0
    if (some.r) {
      k.r <- k[which.r]
      g.r <- g[which.r]
      Ig.r <- Ig[which.r]
      y.r <- y[which.r]
      xtilde.r <- xtilde[which.r]
      onec.r <- onec[which.r]
      cstar.r <- cstar[which.r, , drop = FALSE]
      qg.r <- qg[which.r]
    }
  } else {
    k <- rep(1, n)
    some.r <- FALSE
  }

  # Separate out pools with single Xtilde
  if (errors == "neither") {
    which.i <- NULL
  } else if (errors == "processing") {
    which.i <- which(Ig == 1 & k == 1)
  } else if (errors %in% c("measurement", "both")) {
    which.i <- which(k == 1)
  }
  n.i <- length(which.i)
  some.i <- n.i > 0
  if (some.i) {
    g.i <- g[which.i]
    Ig.i <- Ig[which.i]
    y.i <- y[which.i]
    cstar.i <- cstar[which.i, , drop = FALSE]
    qg.i <- qg[which.i]
    onec.i <- onec[which.i]
    xtilde.i <- unlist(xtilde[which.i])
  }

  # Get indices for parameters being estimated and create labels
  loc.betas <- 1: n.betas
  beta.labels <- paste("beta", c("0", x.varname, c.varnames), sep = "_")

  loc.alphas <- (n.betas + 1): (n.betas + n.alphas)
  alpha.labels <- paste("alpha", c("0", c.varnames), sep = "_")

  loc.b <- n.betas + n.alphas + 1

  if (errors == "neither") {
    theta.labels <- c(beta.labels, alpha.labels, "b")
  } else if (errors == "processing") {
    if (! nondiff_pe) {
      loc.sigsq_p1 <- loc.b + 1
      loc.sigsq_p0 <- loc.b + 2
      theta.labels <- c(beta.labels, alpha.labels,
                        "b", "sigsq_p1", "sigsq_p0")
    } else {
      loc.sigsq_p1 <- loc.sigsq_p0 <- loc.b + 1
      theta.labels <- c(beta.labels, alpha.labels, "b", "sigsq_p")
    }
  } else if (errors == "measurement") {
    if (! nondiff_me) {
      loc.sigsq_m1 <- loc.b + 1
      loc.sigsq_m0 <- loc.b + 2
      theta.labels <- c(beta.labels, alpha.labels,
                        "b", "sigsq_m1", "sigsq_m0")
    } else {
      loc.sigsq_m1 <- loc.sigsq_m0 <- loc.b + 1
      theta.labels <- c(beta.labels, alpha.labels, "b", "sigsq_m")
    }
  } else if (errors == "both") {
    if (! nondiff_pe & ! nondiff_me) {
      loc.sigsq_p1 <- loc.b + 1
      loc.sigsq_p0 <- loc.b + 2
      loc.sigsq_m1 <- loc.b + 3
      loc.sigsq_m0 <- loc.b + 4
      theta.labels <- c(beta.labels, alpha.labels,
                        "b", "sigsq_p1", "sigsq_p0", "sigsq_m1",
                        "sigsq_m0")
    } else if (! nondiff_pe & nondiff_me) {
      loc.sigsq_p1 <- loc.b + 1
      loc.sigsq_p0 <- loc.b + 2
      loc.sigsq_m1 <- loc.sigsq_m0 <- loc.b + 3
      theta.labels <- c(beta.labels, alpha.labels,
                        "b", "sigsq_p1", "sigsq_p0", "sigsq_m")
    } else if (nondiff_pe & ! nondiff_me) {
      loc.sigsq_p1 <- loc.sigsq_p0 <- loc.b + 1
      loc.sigsq_m1 <- loc.b + 2
      loc.sigsq_m0 <- loc.b + 3
      theta.labels <- c(beta.labels, alpha.labels,
                        "b", "sigsq_p", "sigsq_m1", "sigsq_m0")
    } else if (nondiff_pe & nondiff_me) {
      loc.sigsq_p1 <- loc.sigsq_p0 <- loc.b + 1
      loc.sigsq_m1 <- loc.sigsq_m0 <- loc.b + 2
      theta.labels <- c(beta.labels, alpha.labels,
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

    f.alphas <- matrix(f.theta[loc.alphas], ncol = 1)
    f.alpha_0 <- f.alphas[1]
    f.alpha_c <- matrix(f.alphas[-1], ncol = 1)

    f.b <- f.theta[loc.b]

    if (errors == "neither") {
      f.sigsq_p1 <- f.sigsq_p0 <- f.sigsq_m1 <- f.sigsq_m0 <- 0
    }
    if (errors %in% c("processing", "both")) {
      f.sigsq_p1 <- f.theta[loc.sigsq_p1]
      f.sigsq_p0 <- f.theta[loc.sigsq_p0]
    } else {
      f.sigsq_p1 <- f.sigsq_p0 <- 0
    }
    if (errors %in% c("measurement", "both")) {
      f.sigsq_m1 <- f.theta[loc.sigsq_m1]
      f.sigsq_m0 <- f.theta[loc.sigsq_m0]
    } else {
      f.sigsq_m1 <- f.sigsq_m0 <- 0
    }

    if (some.p) {

      # Likelihood for pools with precisely measured X:
      # L = f(Y|X,C) f(X|C_1,...,C_g)

      # P(Y|X,C)
      eta <- gxc.p %*% f.betas + qg.p
      p_y.xc <- (1 + exp(-eta))^(-1)

      # a_i's in X|C_1, ..., C_g ~ Gamma(a_i, b)
      if (some.cs) {
        alphas <- sapply(onec.p, function(x) sum(exp(x %*% f.alphas)))
      } else {
        alphas <- g.p * exp(f.alpha_0)
      }

      # Log-likelihood (non-positive X that generate NAs get excluded!)
      ll.p <- sum(dbinom(x = y.p, size = 1, prob = p_y.xc, log = TRUE) +
                    dgamma(x = x.p, shape = alphas, scale = f.b, log = TRUE))

    } else {
      ll.p <- 0
    }

    # Set skip.rest flag to FALSE
    skip.rest <- FALSE

    if (some.r) {

      # Likelihood for pools with replicates
      # L = \int_X f(Y|X,C) f(Xtilde|X) f(X|C_1,...,C_g) dX

      # Create error vectors
      sigsq_p <- ifelse(y.r > 0, f.sigsq_p1, f.sigsq_p0) * Ig.r
      sigsq_m <- ifelse(y.r > 0, f.sigsq_m1, f.sigsq_m0)

      # a_i's to feed to integral
      if (some.cs) {
        alphas <- sapply(onec.r, function(x) sum(exp(x %*% f.alphas)))
      } else {
        alphas <- g.r * exp(f.alpha_0)
      }

      # Function for integrating out X's
      int.f_i1 <- function(k_i, g_i, Ig_i, y_i, x_i, cstar_i, qg_i, xtilde_i,
                           a_i, sigsq_p_i, sigsq_m_i) {

        # f(Y,Xtilde,X|C_1,...,C_g)
        x_i <- matrix(x_i, nrow = 1)
        f_yxtildex.c <- apply(x_i, 2, function(z) {

          # Transformation
          s_i <- z / (1 - z)

          # P(Y|X,C)
          if (some.cs) {
            p_y.xc <- (1 + exp(-g_i * f.beta_0 - f.beta_x * s_i -
                                 as.vector(t(f.beta_c) %*% cstar_i) - qg_i))^(-1)
          } else {
            p_y.xc <- (1 + exp(-g_i * f.beta_0 - f.beta_x * s_i - qg_i))^(-1)
          }

          # E[log(Xtilde)|X] and V[log(Xtilde|X)]
          Mu_logxtilde.x <- matrix(log(s_i) - 1/2 * (sigsq_p_i + sigsq_m_i), nrow = k_i)
          Sigma_logxtilde.x <- matrix(sigsq_p_i, ncol = k_i, nrow = k_i) +
            diag(rep(sigsq_m_i, k_i))

          # Density
          dbinom(x = y_i, size = 1, prob = p_y.xc) *
            1 / prod(xtilde_i) * dmvnorm(x = log(xtilde_i),
                                         mean = Mu_logxtilde.x,
                                         sigma = Sigma_logxtilde.x) *
            dgamma(x = s_i, shape = a_i, scale = f.b)

        })

        # Back-transformation
        out <- matrix(f_yxtildex.c / (1 - x_i)^2, ncol = ncol(x_i))

      }

      # Get integration tolerance
      if (estimating.hessian) {
        int_tol <- integrate_tol_hessian
      } else if (identical(f.theta, extra.args$start, ignore.environment = TRUE)) {
        int_tol <- integrate_tol_start
      } else {
        int_tol <- integrate_tol
      }

      int.vals <- c()
      for (ii in 1: length(xtilde.r)) {

        # Get values for ith participant
        g_i <- g.r[ii]
        Ig_i <- Ig.r[ii]
        k_i <- k.r[ii]
        y_i <- y.r[ii]
        cstar_i <- cstar.r[ii, ]
        qg_i <- qg.r[ii]
        xtilde_i <- xtilde.r[[ii]]
        a_i <- alphas[ii]
        sigsq_p_i <- sigsq_p[ii]
        sigsq_m_i <- sigsq_m[ii]

        # Try integrating out X with default settings
        int.ii <-
          hcubature(f = int.f_i1, tol = int_tol,
                    lowerLimit = 0, upperLimit = 1,
                    vectorInterface = TRUE,
                    g_i = g_i, Ig_i = Ig_i, k_i = k_i, y_i = y_i,
                    cstar_i = cstar_i, qg_i = qg_i, xtilde_i = xtilde_i,
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

      # Likelihood for pools with single Xtilde:
      # L = \int_X f(Y|X,C) f(Xtilde|X) f(X|C_1,...,C_g) dX

      # Create error vectors
      sigsq_p <- ifelse(y.i > 0, f.sigsq_p1, f.sigsq_p0) * Ig.i
      sigsq_m <- ifelse(y.i > 0, f.sigsq_m1, f.sigsq_m0)

      # a_i's to feed to integral
      if (some.cs) {
        alphas <- sapply(onec.i, function(x) sum(exp(x %*% f.alphas)))
      } else {
        alphas <- g.i * exp(f.alpha_0)
      }

      # Function for integrating out X's
      int.f_i2 <- function(g_i, Ig_i, y_i, x_i, cstar_i, qg_i, xtilde_i, a_i,
                           sigsq_p_i, sigsq_m_i) {

        # Transformation
        s_i <- x_i / (1 - x_i)

        # P(Y|X,C)
        if (some.cs) {
          p_y.xc <- (1 + exp(-g_i * f.beta_0 - f.beta_x * s_i -
                               as.vector(t(f.beta_c) %*% cstar_i) - qg_i))^(-1)
        } else {
          p_y.xc <- (1 + exp(-g_i * f.beta_0 - f.beta_x * s_i - qg_i))^(-1)
        }

        # E[log(Xtilde)|X] and V[log(Xtilde|X)]
        mu_logxtilde.x <- log(s_i) - 1/2 * (sigsq_p_i + sigsq_m_i)
        sigsq_logxtilde.x <- sigsq_p_i + sigsq_m_i

        # f(Y,Xtilde,X|C_1,...,C_g) = f(Y|X,C) f(Xtilde|X) f(X|C_1,...,C_g)
        f_yxtildex.c <- dbinom(x = y_i, size = 1, prob = p_y.xc) *
          1 / xtilde_i * dnorm(x = log(xtilde_i),
                               mean = mu_logxtilde.x,
                               sd = sqrt(sigsq_logxtilde.x)) *
          dgamma(x = s_i, shape = a_i, scale = f.b)

        # Back-transformation
        out <- f_yxtildex.c / (1 - x_i)^2
        return(out)

      }

      # Get integration tolerance
      if (estimating.hessian) {
        int_tol <- integrate_tol_hessian
      } else if (identical(f.theta, extra.args$start, ignore.environment = TRUE)) {
        int_tol <- integrate_tol_start
      } else {
        int_tol <- integrate_tol
      }

      int.vals <- c()
      for (ii in 1: length(xtilde.i)) {

        # Get values for ith participant
        g_i <- g.i[ii]
        Ig_i <- Ig.i[ii]
        y_i <- y.i[ii]
        cstar_i <- cstar.i[ii, ]
        qg_i <- qg.i[ii]
        xtilde_i <- xtilde.i[ii]
        a_i <- alphas[ii]
        sigsq_p_i <- sigsq_p[ii]
        sigsq_m_i <- sigsq_m[ii]

        # Try integrating out X_i with default settings
        int.ii <-
          hcubature(f = int.f_i2, tol = int_tol,
                    lowerLimit = 0, upperLimit = 1,
                    vectorInterface = TRUE,
                    g_i = g_i, Ig_i = Ig_i, y_i = y_i,
                    cstar_i = cstar_i, qg_i = qg_i, xtilde_i = xtilde_i,
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
    if (errors == "neither") {
      extra.args$start <- c(rep(0.01, n.betas + n.alphas), 1)
    } else if (errors == "processing") {
      extra.args$start <- c(rep(0.01, n.betas + n.alphas),
                            rep(1, loc.sigsq_p0 - loc.b + 1))
    } else if (errors %in% c("measurement", "both")) {
      extra.args$start <- c(rep(0.01, n.betas + n.alphas),
                            rep(1, loc.sigsq_m0 - loc.b + 1))
    }
  }
  if (is.null(extra.args$lower)) {
    if (errors == "neither") {
      extra.args$lower <- c(rep(-Inf, n.betas + n.alphas), 1e-3)
    } else if (errors == "processing") {
      extra.args$lower <- c(rep(-Inf, n.betas + n.alphas),
                            rep(1e-3, loc.sigsq_p0 - loc.b + 1))
    } else if (errors %in% c("measurement", "both")) {
      extra.args$lower <- c(rep(-Inf, n.betas + n.alphas),
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
  if (estimate_var) {

    # Estimate Hessian
    hessian.mat <- hessian(f = ll.f, estimating.hessian = TRUE, x0 = theta.hat)
    theta.variance <- try(solve(hessian.mat), silent = TRUE)
    if (class(theta.variance) == "try-error" ||
        ! all(eigen(x = theta.variance, only.values = TRUE)$values > 0)) {

      # Repeatedly divide integrate_tol_hessian by 5 and re-try
      while (integrate_tol_hessian > 1e-15 & fix_posdef) {
        integrate_tol_hessian <- integrate_tol_hessian / 5
        hessian.mat <- hessian(f = ll.f, estimating.hessian = TRUE,
                               x0 = theta.hat)
        theta.variance <- try(solve(hessian.mat), silent = TRUE)
        if (class(theta.variance) != "try-error" &&
            all(eigen(x = theta.variance, only.values = TRUE)$values > 0)) {
          break
        }

      }
    }

    if (class(theta.variance) == "try-error" ||
        ! all(eigen(x = theta.variance, only.values = TRUE)$values > 0)) {

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
