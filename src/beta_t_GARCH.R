
score_beta_t_GARCH <- function(y, lambda_t, nu){
  # return score of lambda of Student-t distribution 
  
  u_t <- ((nu + 1) * y^2) / ((nu - 2) + y^2/lambda_t) - lambda_t
  
  return(u_t)
}


loglik_beta_t_GARCH <- function(y, lambda_t, nu, log = TRUE){
  # compute log-density of Student's t distribution with nu degrees of freedom
  
  ulpdf <- (lgamma((nu+1)/2) - lgamma(nu/2) - (1/2)*log(lambda_t) -
              (1/2)  * log(pi*(nu-2)) - ((nu+1)/2) * log(1 + y^2 / ((nu-2)*lambda_t) ))
  
  if(log != TRUE){
    ulpdf <- exp(ulpdf)
  } 
  
  return(ulpdf)
}



filter_beta_t_GARCH <- function(y, theta){
  # returns the filter of the time-varying scale parameter, the innovation 
  # and the likelihood of the model
  
  T <- length(y)
  
  # Define log likelihoods
  dloglik <- array(data = NA, dim = c(T))
  loglik  <- numeric()
  
  # parameter selections
  omega <- theta[1]
  phi   <- theta[2]
  k     <- theta[3]
  nu    <- theta[4]
  
  # define dynamic scale and score
  lambda_t <- array(data = NA, dim = c(T+1))
  sigma_t  <- array(data = NA, dim = c(T+1))
  u_t      <- array(data = NA, dim = c(T))
  
  # initialize dynamic scale
  lambda_t[1] <- omega/(1-phi)
  sigma_t[1]  <- sqrt(lambda_t[1]/2)
  
  # initialize log density and log likelihood
  dloglik[1] <- loglik_beta_t_GARCH(y[1], lambda_t[1], nu = nu, log = TRUE)
  loglik     <- dloglik[1]
  
  for(t in 2:(T)) {
    # dynamic scale innovations
    u_t[t-1] <- score_beta_t_GARCH(y[t-1], lambda_t[t-1], nu)
    # updating filter                    
    lambda_t[t] <- omega + phi * lambda_t[t-1] + k * u_t[t-1]
    sigma_t[t]  <- sqrt(lambda_t[t]/2)
    
    #Updating log density and log likelihood
    dloglik[t] <- loglik_beta_t_GARCH(y[t], lambda_t = lambda_t[t], nu = nu, log = TRUE)
    loglik     <- loglik + dloglik[t]
    
  }
  
  # output
  lambda_t <- ts(lambda_t, start = start(y), frequency = frequency(y))
  sigma_t  <- ts(sigma_t, start = start(y), frequency = frequency(y))
  u_t      <- ts(u_t, start = start(y), frequency = frequency(y))
  
  out <- list(Dynamic_Scale    = sigma_t,
              Innovation_u_t   = u_t,
              Log_Densities_i  = dloglik,
              Log_Likelihood   = loglik)
  
  return(out)
}


filtered_loglik_beta_t_GARCH <- function(dati, param){
  # the function returns the value of the log-likelihood of the model obtained
  # after filtering the time-varying scale parameter with some prespecified
  # static parameters
  
  # parameter elections 
  omega <- param[1]
  phi   <- param[2]
  k     <- param[3]
  nu    <- param[4]
  
  theta_new <- c(omega, phi, k, nu)
  
  # extract log-likelihood
  fitness <- filter_beta_t_GARCH(dati, theta_new)$Log_Likelihood
  
  if(is.na(fitness) | !is.finite(fitness)) fitness <- -1e10
  if(fitness != fitness) fitness <- -1e10
  
  return(-fitness)
} 


estimator_beta_t_GARCH <- function(dati, param){
  # computes the optimal static parameters iteratively, starting with a set of
  # given values and by maximizing the log-likelihood obtained by filtering the 
  # time-varying scale parameter
  
  # parameter elections 
  omega <- param[1]
  phi   <- param[2]
  k     <- param[3]
  nu    <- param[4]
  
  theta_st <- c(omega, phi, k, nu)
  
  # parameter constraints
  lower <- c(-Inf, -0.999, -2, 2.099)
  upper <- c(Inf,  0.999,  2, 300)
  
  # Optimize every Filters w/nlminb 
  optimizer <- suppressWarnings(nlminb(start = theta_st, objective = filtered_loglik_beta_t_GARCH, 
                                       dati  = dati, gradient = NULL, 
                                       control = list(trace = 0), hessian = NULL,
                                       lower = lower, upper = upper))
  
  # Save the optimized parameters Dynamic Location
  omega_opt <- optimizer$par[1]  
  phi_opt   <- optimizer$par[2]
  k_opt     <- optimizer$par[3]
  nu_opt    <- optimizer$par[4]
  
  # Create a vector with all the optimized parameters
  theta_opt <- c(omega_opt, phi_opt, k_opt, nu_opt)
  
  # Create a list with all the optimized parameters
  theta_list <- list(omega = omega_opt,
                     phi   = phi_opt,
                     k     = k_opt,
                     nu    = nu_opt)
  
  # Output
  out <- list(theta_list = theta_list,
              theta      = theta_opt,
              optimizer  = optimizer)
  
  return(out) 
}

