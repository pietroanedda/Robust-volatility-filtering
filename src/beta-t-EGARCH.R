
score_beta_t_EGARCH <- function(y, lambda_t, nu){
  # return score of variance sigma^2 of Student-t distribution 
  
  u_t <- ((nu + 1) * y^2) / ((nu - 2) * exp(lambda_t) + y^2) - 1
  
  return(u_t)
}


loglik_beta_t_EGARCH <- function(y, lambda_t, nu, log = TRUE){
  # compute Student-t log-density
  
  ulpdf <- (lgamma((nu+1)/2) - lgamma(nu/2) - (1/2) * lambda_t -
              (1/2)  * log(pi*(nu-2)) - ((nu+1)/2) * log(1 + y^2 / ((nu-2)*exp(lambda_t)) ))
  
  if(log != TRUE){
    ulpdf <- exp(ulpdf)
  } 
  
  return(ulpdf)
}



filter_beta_t_EGARCH <- function(y, theta){
  # returns the filter of the time-varying scale parameter 
  
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
  sigma_t <- array(data = NA, dim = c(T+1))
  u_t      <- array(data = NA, dim = c(T))
  
  # initialize dynamic scale
  lambda_t[1] <- omega/(1-phi)
  sigma_t[1] <- exp(lambda_t[1]/2)
  
  # initialize log likelihood
  dloglik[1] <- loglik_beta_t_EGARCH(y[1], lambda_t[1], nu = nu, log = TRUE)
  loglik     <- dloglik[1]
  
  for(t in 2:(T)) {
    # dynamic scale innovations
    u_t[t-1] <- score_beta_t_EGARCH(y[t-1], lambda_t[t-1], nu)
    # updating filter                    
    lambda_t[t] <- omega + phi * lambda_t[t-1] + k * u_t[t-1]
    sigma_t[t] <- exp(lambda_t[t]/2)
    
    if(t < (T+1)){
      #Updating likelihoods
      dloglik[t] <- loglik_beta_t_EGARCH(y[t], lambda_t = lambda_t[t], nu = nu, log = TRUE)
      loglik     <- loglik + dloglik[t]
    }
  }
  
  # output
  lambda_t <- ts(lambda_t, start = start(y), frequency = frequency(y))
  sigma_t <- ts(sigma_t, start = start(y), frequency = frequency(y))
  u_t <- ts(u_t, start = start(y), frequency = frequency(y))
  
  out <- list(Dynamic_Scale    = sigma_t,
              Innovation_u_t   = u_t,
              Log_Densities_i  = dloglik,
              Log_Likelihood   = loglik)
  
  return(out)
}


filtered_loglik_beta_t_EGARCH <- function(dati, param){
  
  # Parameter Selections Dynamic Location
  omega <- param[1]
  phi   <- param[2]
  k     <- param[3]
  nu    <- param[4]
  
  # Create a new vector with the parameters
  theta_new <- c(omega, phi, k, nu)
  
  # Fitness Functions
  fitness <- filter_beta_t_EGARCH(dati, theta_new)$Log_Likelihood
  
  if(is.na(fitness) | !is.finite(fitness)) fitness <- -1e10
  if(fitness != fitness) fitness <- -1e10
  
  return(-fitness)
} 


estimator_beta_t_EGARCH <- function(dati, param){
  
  # Parameter Selections Dynamic Location
  omega <- param[1]
  phi   <- param[2]
  k     <- param[3]
  nu    <- param[4]
  
  # Create a vector with the parameters
  theta_st <- c(omega, phi, k, nu)
  
  # Take Bounds
  lower <- c(-Inf, -0.999, -2, 2.099)
  upper <- c(Inf,  0.999,  2, 300)
 
  # Optimize every Filters w/nlminb 
  optimizer <- suppressWarnings(nlminb(start = theta_st, objective = filtered_loglik_beta_t_EGARCH, 
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

