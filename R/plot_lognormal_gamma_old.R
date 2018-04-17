#' Generate Data for Lognormal or Gamma Regression Model and Plot Fitted Models
#'
#' Useful for assessing whether lognormal and Gamma models are similar for
#' various parameter values.
#'
#' If \code{b} is not specified, data are generated:
#'
#' C ~ N(0, 1)
#' X|C ~ Lognormal(beta_0 + beta_c C, sigsq)
#'
#' If \code{b} is specified, data are generated:
#'
#' C ~ N(0, 1)
#' X|C ~ Gamma(exp(beta_0 + beta_c C), b)
#'
#' Both models are fit using maximum likelihood, and a plot is generated
#' showing E(X) vs. C for each.
#'
#' @param n Numeric value.
#' @param beta_0 Numeric value.
#' @param beta_c Numeric value.
#' @param sigsq Numeric value.
#' @param b Numeric value
#'
#'
#' @return Plot of E(X) vs. X according to each fitted model.
#'
#'
#'@export
plot_lognormal_gamma <- function(n = 1000,
                                 beta_0 = 0,
                                 beta_c = 0.5,
                                 sigsq = 0.1,
                                 b = NULL) {

  # Generate data
  c <- rnorm(n = n)
  if (is.null(b)) {
    x <- rlnorm(n = n, meanlog = beta_0 + beta_c * c, sdlog = sqrt(sigsq))
  } else {
    x <- rgamma(n = n, shape = exp(beta_0 + beta_c * c), scale = b)
  }

  # Fit models
  fit.gamma1 <- glm(x ~ c, family = Gamma("identity"))
  fit.gamma2 <- glm(x ~ c, family = Gamma(""))


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
    logOR <- (sigsq_1 - sigsq_0) / (sigsq_1 * sigsq_0) +
      x * (sigsq_1 - sigsq_0) / (2 * sigsq_1 * sigsq_0) +
      (sigsq_0 * (gamma_0 + gamma_y) - sigsq_1 * gamma_0) /
      (sigsq_1 * sigsq_0)
    df <- data.frame(x = x, logOR = logOR)

    # Calculate confidence bands
    if (! is.null(varcov)) {

      ses <- sapply(x, function(x) {
        fprime <- matrix(c((sigsq_0 - sigsq_1) / (sigsq_1 * sigsq_0),
                           sigsq_0 / (sigsq_1 * sigsq_0),
                           1 / sigsq_1^2 + x / (2 * sigsq_1^2) -
                             (gamma_0 + gamma_y) / sigsq_1^2,
                           -1 / sigsq_0^2 - x / (2 * sigsq_0^2) +
                             gamma_0 / sigsq_0^2),
                         nrow = 1)
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
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
        ylim(min(df$lower), max(df$upper))

    }

  } else if (is.numeric(cvals)) {

    # 1 set of covariate values - plot curve and confidence bands (if possible)

    # Calculate log-OR's
    logOR <- (sigsq_1 - sigsq_0) / (sigsq_1 * sigsq_0) +
      x * (sigsq_1 - sigsq_0) / (2 * sigsq_1 * sigsq_0) +
      (sigsq_0 * (gamma_0 + sum(gamma_c * cvals) + gamma_y) -
         sigsq_1 * (gamma_0 + sum(gamma_c * cvals))) /
      (sigsq_1 * sigsq_0)
    df <- data.frame(x = x, logOR = logOR)

    # Calculate confidence bands
    if (! is.null(varcov)) {

      ses <- sapply(x, function(x) {
        fprime <- matrix(c((sigsq_0 - sigsq_1) / (sigsq_1 * sigsq_0),
                           sigsq_0 / (sigsq_1 * sigsq_0),
                           cvals * (sigsq_0 - sigsq_1) / (sigsq_1 * sigsq_0),
                           1 / sigsq_1^2 + x / (2 * sigsq_1^2) -
                             (gamma_0 + sum(gamma_c * cvals) + gamma_y) / sigsq_1^2,
                           -1 / sigsq_0^2 - x / (2 * sigsq_0^2) +
                             (gamma_0 + sum(gamma_c * cvals)) / sigsq_0^2),
                         nrow = 1)
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
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
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
      logOR <- (sigsq_1 - sigsq_0) / (sigsq_1 * sigsq_0) +
        x * (sigsq_1 - sigsq_0) / (2 * sigsq_1 * sigsq_0) +
        (sigsq_0 * (gamma_0 + sum(gamma_c * cvals.ii) + gamma_y) -
           sigsq_1 * (gamma_0 + sum(gamma_c * cvals.ii))) /
        (sigsq_1 * sigsq_0)
      df <- bind_rows(df, data.frame(Covariates = ii, x = x, logOR = logOR))

      # Calculate confidence bands
      if (! is.null(varcov) & set_panels) {

        ses <- sapply(x, function(x) {
          fprime <- matrix(c((sigsq_0 - sigsq_1) / (sigsq_1 * sigsq_0),
                             sigsq_0 / (sigsq_1 * sigsq_0),
                             cvals.ii * (sigsq_0 - sigsq_1) / (sigsq_1 * sigsq_0),
                             1 / sigsq_1^2 + x / (2 * sigsq_1^2) -
                               (gamma_0 + sum(gamma_c * cvals.ii) + gamma_y) / sigsq_1^2,
                             -1 / sigsq_0^2 - x / (2 * sigsq_0^2) +
                               (gamma_0 + sum(gamma_c * cvals.ii)) / sigsq_0^2),
                           nrow = 1)
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
        facet_grid(facets = . ~ Covariates) +
        geom_line() +
        geom_hline(yintercept = 0, linetype = 2) +
        labs(title = paste("Log-OR vs.", xname),
             y = "Log-OR",
             x = xname) +
        ylim(min(logOR), max(logOR)) +
        theme_bw()

      if (! is.null(varcov)) {

        p <- p +
          geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
          ylim(min(df$lower), max(df$upper))

      }

    } else {

      p <- ggplot(df, aes(x, logOR, group = Covariates, color = Covariates)) +
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
