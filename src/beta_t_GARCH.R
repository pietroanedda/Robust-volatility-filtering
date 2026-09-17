# FUNCTIONS TO FIT A BETA-T-GARCH MODEL

# the procedure starts by computing a general expression for the loglikelihood of the model
# the loglikelihood is then maximized via a quasi-Newton algorithm, obtaining optimal parameters
# given the optimal parameters, the TVP is filtered and the likelihood, the conditional score
# and the standardized residuals are computed


score_beta_t_GARCH <- function(y, lambda, nu){
  # return score of lambda of Student-t distribution 
  # lambda is the vector of the time-varying parameter
  
  u_t <- ((nu + 1) * y^2) / ((nu - 2) + y^2/lambda) - lambda
  return(u_t)
}


loglik_beta_t_GARCH <- function(y, lambda, nu, log = TRUE){
  # compute log-density of Student's t distribution with nu degrees of freedom
  # by default returns log-likelihood, if likelihood is required set log=FALSE
  
  ulpdf <- (lgamma((nu+1)/2) - lgamma(nu/2) - (1/2)*log(lambda) -
              (1/2)  * log(pi*(nu-2)) - ((nu+1)/2) * log(1 + y^2 / ((nu-2)*lambda) ))
  
  # if log-likelihood not required, compute likelihood
  if(log != TRUE){
    ulpdf <- exp(ulpdf)
  } 
  
  return(ulpdf)
}


filter_beta_t_GARCH <- function(y, theta){
  # returns the filter of the time-varying scale parameter, the innovation 
  # series and the likelihood of the model
  # the vector theta must contain: intercept, autoregressive coeff, score coeff 
  # and degrees of freedom, in that order
  
  T <- length(y)
  
  # Define log likelihoods
  dloglik <- numeric(T)
  loglik  <- 0
  
  # parameter selections from theta vector
  omega <- theta[1]
  phi   <- theta[2]
  k     <- theta[3]
  nu    <- theta[4]
  
  # define dynamic scale, volatility and conditional score
  lambda <- numeric(T)
  sigma  <- numeric(T)
  u      <- numeric(T)
  
  # initialize dynamic scale
  # take unconditional expectation of the autoregressive process
  lambda[1] <- omega/(1-phi)
  sigma[1]  <- sqrt(lambda[1])
  
  # initialize log density and log likelihood
  dloglik[1] <- loglik_beta_t_GARCH(y[1], lambda[1], nu = nu, log = TRUE)
  loglik     <- dloglik[1]
  
  for(t in 2:(T)) {
    # dynamic scale innovations
    u[t-1] <- score_beta_t_GARCH(y[t-1], lambda[t-1], nu)
    
    # updating filter                    
    lambda[t] <- omega + phi * lambda[t-1] + k * u[t-1]
    sigma[t]  <- sqrt(lambda[t])
    
    #Updating log density series and log likelihood
    dloglik[t] <- loglik_beta_t_GARCH(y[t], lambda=lambda[t], nu=nu, log=TRUE)
    loglik     <- loglik + dloglik[t]
  }
  
  # compute conditional score for last period
  u[T] <- score_beta_t_GARCH(y[T], lambda[T], nu)
  
  # output
  lambda <- ts(lambda, start = start(y), frequency = frequency(y))
  sigma  <- ts(sigma, start = start(y), frequency = frequency(y))
  u      <- ts(u, start = start(y), frequency = frequency(y))
  
  out <- list(Dynamic_Scale    = sigma,
              Innovation_u_t   = u,
              Log_Densities_i  = dloglik,
              Log_Likelihood   = loglik)
  
  return(out)
}


filtered_loglik_beta_t_GARCH <- function(param, data){
  # the function returns the value of the log-likelihood of the model obtained
  # after filtering the time-varying scale parameter with some prespecified
  # static parameters
  
  # extract log-likelihood
  fitness <- filter_beta_t_GARCH(data, param)$Log_Likelihood
  
  if(is.na(fitness) | !is.finite(fitness)) fitness <- -1e10
  if(fitness != fitness) fitness <- -1e10
  
  return(-fitness)
} 


estimator_beta_t_GARCH <- function(data, param){
  # computes the optimal static parameters iteratively, starting with a set of
  # given values and by maximizing the log-likelihood obtained by filtering the 
  # time-varying scale parameter
  
  # parameter constraints for Beta-t-GARCH (existence and stationarity)
  lower <- c(1e-8, 0, 0, 2.099)
  upper <- c(Inf,  0.99,  10, 100)
  
  # Optimize parameters w/nlminb (quasi-Newton, similar to L-BFGS-B)
  optimizer <- suppressWarnings(nlminb(start = param, objective = filtered_loglik_beta_t_GARCH, 
                                       data  = data, gradient = NULL, 
                                       control = list(trace = 0), hessian = NULL,
                                       lower = lower, upper = upper))
  
  # Create a list with all the optimized parameters
  theta_list <- list(omega = optimizer$par[1],
                     phi   = optimizer$par[2],
                     k     = optimizer$par[3],
                     nu    = optimizer$par[4])
  
  # Output 
  out <- list(theta_list = theta_list,
              theta      = optimizer$par,
              optimizer  = optimizer)
  
  return(out) 
}

