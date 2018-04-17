#' Fit Lognormal Regression for Y vs. Covariates
#' 
#' Uses maximum likelihood to fit 
#' Y|\strong{X} ~ Lognormal(beta_0 + \strong{beta_x}^T \strong{X}, sigsq)
#'
#' @param y Numeric vector.
#' @param x Numeric vector or matrix. If \code{NULL}, model reduces to marginal
#' lognormal model Y ~ Lognormal(exp(beta_0), sigsq).
#' @param var Logical value for whether to return Hessian-based
#' variance-covariance matrix.
#'
#' @return List of parameter estimates, variance-covariance matrix (if
#' requested), AIC, and \code{\link[stats]{nlminb}} object.
#'
#'
#' @examples
#' # Generate data
#' set.seed(123)
#' x <- rnorm(1000)
#' y <- rlnorm(1000, meanlog = 0.5 + 0.25 * x, sdlog = 0.5)
#'
#' # Fit model
#' fit <- lognormal(y = y, x = x)
#' fit$theta.hat
#' fit$varcov
#' fit$aic
#'
#' # Plot E(Y) vs. X according to model fit
#' plot(x, y, main = "Lognormal Model for Y vs. X")
#' xvals <- seq(min(x), max(x), 0.01)
#' yvals <- exp(fit$theta.hat[1] + fit$theta.hat[2] * xvals + fit$theta.hat[3] / 2)
#' points(xvals, yvals, type = "l")
#'
#'
#'@export
lognormal <- function(y, x = NULL, var = TRUE) {
  
  # Design matrix
  n <- length(y)
  onex <- cbind(rep(1, n), x)
  p <- ncol(onex)
  
  # Labels
  labs <- c(paste("x", 0: (p - 1), sep = ""), "b")
  
  # Log-likelihood function
  ll.f <- function(f.theta) {
    
    f.betas <- f.theta[1: p]
    f.sigsq <- f.theta[p + 1]
    ll <- sum(dlnorm(x = y, log = TRUE,
                     meanlog = onex %*% f.betas, 
                     sdlog = sqrt(f.sigsq)))
    return(-ll)
    
  }
  
  # Maximize log-likelihood
  llmax <- nlminb(objective = ll.f,
                  start = c(rep(0, p), 1),
                  lower = c(rep(-Inf, p), 1e-6))
  theta.hat <- llmax$par
  names(theta.hat) <- labs
  ret.list <- list(theta.hat = theta.hat)
  
  # Estimate variance-covariance matrix
  if (var) {
    hessian.mat <- pracma::hessian(f = ll.f, x0 = theta.hat)
    varcov <- solve(hessian.mat)
    colnames(varcov) <- rownames(varcov) <- labs
    ret.list$varcov <- varcov
  }
  
  # Add AIC and nlminb object
  ret.list$aic <- 2 * (p + 1 + llmax$objective)
  ret.list$nlminb.object <- llmax
  
  # Return ret.list
  return(ret.list)
  
}
