#' Plot Log-OR vs. X for Normal Discriminant Function Approach
#'
#' When \code{\link{p_ndfa}} is fit with \code{constant_or = FALSE}, the
#' log-OR for X depends on the value of X (and covariates, if any). This
#' function plots the log-OR vs. X for one or several sets of covariate values.
#'
#'
#' @param estimates Numeric vector of point estimates for
#' \code{(gamma_0, gamma_y, gamma_c^T, sigsq)}.
#' @param varcov Numeric matrix with variance-covariance matrix for
#' \code{estimates}. If \code{NULL}, 95\% confidence bands are omitted.
#' @param xrange Numeric vector specifying range of X values to plot.
#' @param xname Character vector specifying name of X variable, for
#' plot title and x-axis label.
#' @param cvals Numeric vector or list of numeric vectors specifying covariate
#' values to use in log-odds ratio calculations.
#' @param set_labels Character vector of labels for the sets of covariate
#' values. Only used if \code{cvals} is a list.
#' @param set_panels Logical value for whether to use separate panels for each
#' set of covariate values, as opposed to using different colors on a single
#' plot.
#'
#'
#' @return Plot of log-OR vs. X generated by \code{\link[ggplot2]{ggplot}}.
#'
#'
#' @examples
#' # Fit discriminant function model for poolwise X vs. (Y, C), without assuming
#' # a constant log-OR. Note that data were generated with a constant log-OR of
#' # 0.5.
#' data(dat_p_ndfa)
#' dat <- dat_p_ndfa$dat
#' fit <- p_ndfa(
#'   g = dat$g,
#'   y = dat$numcases,
#'   xtilde = dat$x,
#'   c = dat$c,
#'   errors = "neither",
#'   constant_or = FALSE
#' )
#'
#' # Plot estimated log-OR vs. X, holding C fixed at the sample mean.
#' p <- plot_ndfa(
#'   estimates = fit$estimates,
#'   varcov = fit$theta.var,
#'   xrange = range(dat$x[dat$g == 1]),
#'   cvals = mean(dat$c / dat$g)
#' )
#' p
#'
#'
#'@export
plot_ndfa <- function(estimates,
                      varcov = NULL,
                      xrange,
                      xname = "X",
                      cvals = NULL,
                      set_labels = NULL,
                      set_panels = TRUE) {

  # Extract parameter estimates
  names_estimates <- names(estimates)

  loc.gammas <- which(substr(names_estimates, 1, 6) == "gamma_")
  loc.sigsq_1 <- which(names_estimates == "sigsq_1")
  loc.sigsq_0 <- which(names_estimates == "sigsq_0")

  gammas <- estimates[loc.gammas]
  gamma_0 <- gammas[1]
  gamma_y <- gammas[2]
  gamma_c <- gammas[-c(1, 2)]
  sigsq_1 <- estimates[loc.sigsq_1]
  sigsq_0 <- estimates[loc.sigsq_0]

  # Subset useful part of variance-covariance matrix
  locs <- c(loc.gammas, loc.sigsq_1, loc.sigsq_0)
  varcov <- varcov[locs, locs]

  # Create X vector
  x <- seq(xrange[1], xrange[2], (xrange[2] - xrange[1]) / 500)

  if (is.null(cvals)) {

    # No-covariate case - plot curve and confidence bands (if possible)

    # Calculate log-OR's
    logOR <- gamma_y / sigsq_1 +
      (1 / sigsq_0 - 1 / sigsq_1) *
      (x - gamma_0 + 1 / 2)
    df <- data.frame(x = x, logOR = logOR)

    # Calculate confidence bands
    if (! is.null(varcov)) {

      ses <- sapply(x, function(x) {

        fprime <- matrix(c(
          1 / sigsq_1 - 1 / sigsq_0,
          1 / sigsq_1,
          1 / sigsq_1^2 * (x - gamma_0 - gamma_y + 1/2),
          -1 / sigsq_0^2 * (x - gamma_0 + 1/2)
        ), nrow = 1)
        sqrt(fprime %*% varcov %*% t(fprime))

      })
      df$lower <- logOR - qnorm(0.975) * ses
      df$upper <- logOR + qnorm(0.975) * ses

    }

    # Create plot
    p <- ggplot(df, aes(x, logOR)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(title = paste("Estimated Log-OR vs.", xname),
           y = "Log-OR",
           x = xname) +
      ylim(min(logOR), max(logOR)) +
      theme_bw()

    # Add confidence bands
    if (! is.null(varcov)) {

      p <- p +
        geom_ribbon(aes_string(ymin = "lower", ymax = "upper"),
                    alpha = 0.2) +
        ylim(min(df$lower), max(df$upper))

    }

  } else if (is.numeric(cvals)) {

    # 1 set of covariate values - plot curve and confidence bands (if possible)

    # Calculate log-OR's
    logOR <- gamma_y / sigsq_1 +
      (1 / sigsq_0 - 1 / sigsq_1) *
      (x - gamma_0 - sum(gamma_c * cvals) + 1 / 2)
    df <- data.frame(x = x, logOR = logOR)

    # Calculate confidence bands
    if (! is.null(varcov)) {

      ses <- sapply(x, function(x) {

        fprime <- matrix(c(
          1 / sigsq_1 - 1 / sigsq_0,
          1 / sigsq_1,
          (1 / sigsq_1 - 1 / sigsq_0) * cvals,
          1 / sigsq_1^2 * (x - gamma_0 - gamma_y - sum(gamma_c * cvals) + 1/2),
          -1 / sigsq_0^2 * (x - gamma_0 - sum(gamma_c * cvals) + 1/2)
        ), nrow = 1)
        sqrt(fprime %*% varcov %*% t(fprime))

      })

      df$lower <- logOR - qnorm(0.975) * ses
      df$upper <- logOR + qnorm(0.975) * ses

    }

    # Create plot
    p <- ggplot(df, aes(x, logOR)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(title = paste("Log-OR vs.", xname),
           y = "Log-OR",
           x = xname) +
      ylim(min(logOR), max(logOR)) +
      theme_bw()

    # Add confidence bands
    if (! is.null(varcov)) {

      p <- p +
        geom_ribbon(aes_string(ymin = "lower", ymax = "upper"),
                    alpha = 0.2) +
        ylim(min(df$lower), max(df$upper))

    }

  } else if (is.list(cvals)) {

    # Multiple sets of covariate values

    # Create labels for covariate sets
    if (is.null(set_labels)) {
      cnames <- substr(names(gamma_c), start = 7, stop = 100)
      set_labels <- sapply(cvals, function(x) paste(cnames, "=", x, collapse = ", "))
    }

    # Loop through covariate sets and calculate log-OR's for each
    df <- NULL
    for (ii in 1: length(cvals)) {

      # Calculate log-OR's
      cvals.ii <- cvals[[ii]]
      logOR <- gamma_y / sigsq_1 +
        (1 / sigsq_0 - 1 / sigsq_1) *
        (x - gamma_0 - sum(gamma_c * cvals.ii) + 1 / 2)

      df <- dplyr::bind_rows(df, data.frame(Covariates = ii, x = x, logOR = logOR))

      # Calculate confidence bands
      if (! is.null(varcov) & set_panels) {

        ses <- sapply(x, function(x) {
          fprime <- matrix(c(
            1 / sigsq_1 - 1 / sigsq_0,
            1 / sigsq_1,
            (1 / sigsq_1 - 1 / sigsq_0) * cvals.ii,
            1 / sigsq_1^2 * (x - gamma_0 - gamma_y - sum(gamma_c * cvals.ii) + 1/2),
            -1 / sigsq_0^2 * (x - gamma_0 - sum(gamma_c * cvals.ii) + 1/2)
          ), nrow = 1)
          sqrt(fprime %*% varcov %*% t(fprime))

        })

        df$lower <- logOR - qnorm(0.975) * ses
        df$upper <- logOR + qnorm(0.975) * ses

      }

    }
    df$Covariates <- factor(df$Covariates, levels = 1: length(cvals), labels = set_labels)

    # Create plot
    if (set_panels) {

      p <- ggplot(df, aes(x, logOR)) +
        facet_grid(reformulate("Covariates", ".")) +
        geom_line() +
        geom_hline(yintercept = 0, linetype = 2) +
        labs(title = paste("Log-OR vs.", xname),
             y = "Log-OR",
             x = xname) +
        ylim(min(logOR), max(logOR)) +
        theme_bw()

      if (! is.null(varcov)) {

        p <- p +
          geom_ribbon(aes_string(ymin = "lower", ymax = "upper"),
                      alpha = 0.2) +
          ylim(min(df$lower), max(df$upper))

      }

    } else {

      p <- ggplot(df, aes_string(x = "x",
                                 y = "logOR",
                                 group = "Covariates",
                                 color = "Covariates")) +
        geom_line() +
        geom_hline(yintercept = 0, linetype = 2) +
        labs(title = paste("Log-OR vs.", xname),
             y = "Log-OR",
             x = xname) +
        ylim(min(logOR), max(logOR)) +
        theme_bw()

    }

  }

  # Plot
  p

}
