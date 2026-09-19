acf_test <- function(y, lags = 10, type="Ljung-Box") {
  # computes the p-values of a test on the autocorrelation of a given time series up to a given lag
  # default test is Ljung-Box test; for each lag, the p-values are plotted
  
  types <- c("Box-Pierce", "Ljung-Box")
  
  stopifnot(
    "'y' must be a numeric vector or a time series" = (is.numeric(y) && is.null(dim(y))) || is.ts(y),
    "'lags' must be a positive integer" = is.numeric(lags) && length(lags) == 1 && lags %% 1 == 0 && lags > 0,
    "'type' must be 'Ljung-Box' (default) or 'Box-Pierce'." = type %in% types
  )
  
  pvalues <- sapply(1:lags, function(k) {
    p  <- Box.test(y, lag=k, type=type)$p.value
    return(p)
  })
  
  plot(1:lags, pvalues, type="p", ylim=c(0, 1), xlab="Lag", ylab="P-value") 
  abline(h = 0.05, lty = 2, col="blue")
}


qqplot_t <- function(y, dof) {
  # plots the qq-plot for a series vs a theoretical Student's t distribution
  # specify the degrees of freedom of the theoretical distribution
  
  stopifnot(
    "'y' must be a numeric vector or a time series" = (is.numeric(y) && is.null(dim(y))) || is.ts(y),
    "'dof' must be a positive integer" = is.numeric(dof) && length(dof) == 1 && dof %% 1 == 0 && dof > 0
  )
  
  # computes the theoretical quantiles 
  t_quantiles <- qt(ppoints(length(y)), df = dof)
  
  qqplot(t_quantiles, y, main="Student's t Q-Q Plot", xlab='Theoretical Quantiles', ylab='Sample Quantiles')
  qqline(y, distribution = function(p) qt(p, df = dof), col = "blue", lwd = 2, lty = 2)
}


IC <- function(lik, param, n) {
  # computes the AIC and BIC of a model; requires the log-likelihood, 
  # a vector of parameters and the number of observations as inputs
  
  stopifnot(
    "'lik' must be a scalar number" = is.numeric(lik) && length(lik) == 1,
    "'param' must be a numeric vector" = is.numeric(param) && is.null(dim(param)),
    "'n' must be a positive integer" = is.numeric(n) && length(n) == 1 && n %% 1 == 0 && n > 0 
  )
  
  aic <- 2 * length(param) - 2*lik
  bic <- log(n) * length(param) - 2*lik
  
  out <- list(AIC=aic, BIC=bic)
  return(out)
}
