# FUNCTIONS TO FIT A GAUSSIAN GARCH MODEL

# the procedure starts by computing a general expression for the loglikelihood of the model
# the loglikelihood is then maximized via a quasi-Newton algorithm, obtaining optimal parameters
# given the optimal parameters, the TVP is filtered and the likelihood, the conditional score and the standardized residuals are computed


loglik_GARCH <- function(y, lambda, log = TRUE){
  # compute log-density of normal distribution 
  # by default returns log-likelihood, if likelihood is required set log=FALSE
  
  ulpdf <- - (1/2)*log(2*pi) - log(sqrt(lambda)) - (1/(2 * lambda)) * y^2 
  
  # if log-likelihood not required, compute likelihood
  if(log != TRUE){
    ulpdf <- exp(ulpdf)
  } 
  
  return(ulpdf)
}


filter_GARCH <- function(y, theta){
  # returns the filter of the time-varying scale parameter, the innovation series and the likelihood of the model
  # the vector theta must contain: omega, alpha, beta, in that order
  # the recursion is specified as lambda_{t+1} = omega + alpha * y_t^2 + beta * lambda_t
  
  n <- length(y)
  
  # Define log likelihoods
  dloglik <- numeric(n)
  loglik  <- 0
  
  # parameter selections from theta vector
  omega <- theta[1]
  alpha <- theta[2]
  beta  <- theta[3]
  
  # define dynamic scale, volatility and MDS
  lambda  <- numeric(n)
  sigma   <- numeric(n)
  v       <- numeric(n)
  std_res <- numeric(n)
  
  # initialize dynamic scale, MDS and standardized residuals
  # take unconditional expectation of the autoregressive process
  lambda[1]  <- omega/(1 - alpha - beta)
  sigma[1]   <- sqrt(lambda[1])
  v[1]       <- (y[1]^2) - lambda[1]
  std_res[1] <- y[1]/sigma[1]
  
  # initialize log density and log likelihood
  dloglik[1] <- loglik_GARCH(y[1], lambda[1], log = TRUE)
  loglik     <- dloglik[1]
  
  for(t in 2:(n)) {
    # updating filter                    
    lambda[t] <- omega + beta * lambda[t-1] + alpha * (y[t-1]^2)
    sigma[t]  <- sqrt(lambda[t])
    std_res[t] <- y[t]/sigma[t]
    
    # martingale difference sequence
    v[t] <- y[t]^2 - lambda[t]
    
    #Updating log density series and log likelihood
    dloglik[t] <- loglik_GARCH(y[t], lambda=lambda[t], log=TRUE)
    loglik     <- loglik + dloglik[t]
  }
  
  # output
  lambda  <- ts(lambda, start = start(y), frequency = frequency(y))
  sigma   <- ts(sigma, start = start(y), frequency = frequency(y))
  v       <- ts(v, start = start(y), frequency = frequency(y))
  std_res <- ts(std_res, start = start(y), frequency = frequency(y))
  
  out <- list(Dynamic_Scale  = sigma,
              Innovation_u_t = v,
              Log_Likelihood = loglik,
              std_residuals  = std_res)
  
  return(out)
}


# the nlminb function takes only box-constraints as inputs, not linear inequality constraints
# the stationarity constraint alpha + beta < 1 cannot be imposed in the optimization
# this problem can be bypassed by adding a penalty to the log-likelihood if the constraint is not observed
# hence the static parameters are estimated through a penalized maximum likelihood approach


conditional_loglik_GARCH <- function(param, data){
  # the function returns the value of the log-likelihood of the model obtained after 
  #filtering the time-varying scale parameter with some prespecified static parameters
  
  # adding penalty to ensure stationarity
  if ((param[2] + param[3]) >= 1) {
    return(1e10)
  }
  
  # extract conditional log-likelihood
  fitness <- filter_GARCH(data, param)$Log_Likelihood
  
  if(is.na(fitness) | !is.finite(fitness)) fitness <- -1e10
  if(fitness != fitness) fitness <- -1e10
  
  return(-fitness)
} 


estimator_GARCH <- function(data, param){
  # computes the optimal static parameters iteratively, starting with a set of
  # given values and by maximizing the log-likelihood obtained by filtering the 
  # time-varying scale parameter
  
  # parameter constraints for GARCH 
  lower <- c(1e-8, 0, 0)
  upper <- c(Inf,  0.999,  0.999)
  
  # Optimize parameters w/nlminb (quasi-Newton, similar to L-BFGS-B)
  optimizer <- nlminb(start = param, objective = conditional_loglik_GARCH, 
                      data  = data, gradient = NULL, 
                      control = list(trace = 0), hessian = NULL,
                      lower = lower, upper = upper)
  
  # check for algorithm convergence
  if(optimizer$convergence != 0) warning("nlminb: ", optimizer$message)
  
  # Create a list with all the optimized parameters
  theta_list <- list(omega = optimizer$par[1],
                     alpha = optimizer$par[2],
                     beta  = optimizer$par[3])
  
  # Output 
  out <- list(theta_list = theta_list,
              theta      = optimizer$par,
              optimizer  = optimizer)
  
  return(out) 
}
