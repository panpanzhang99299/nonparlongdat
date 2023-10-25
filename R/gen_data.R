gen_data <- function(X.timevar = FALSE, aux.timevar = FALSE,
                     X.cont = TRUE, aux.cont = TRUE){
  
  if ((X.cont == TRUE) && (aux.cont == TRUE)){
    
    ## Auxiliary variables
    mean_X_aux <- 3
    sd_X_aux <- 0.1
    
    if (aux.timevar == FALSE){
      
      if (X.timevar == FALSE){
        
        ## X: continuous, time-invariant; aux: continuous, time-invariant
        
        X_aux <- rnorm(n, mean = mean_X_aux, sd = sd_X_aux)
        X_aux_vec <- rep(X_aux, m)
        x <- rnorm(n, mean = 0, sd = 1)
        x.perp <- residuals(lm(x ~ X_aux))
        X <- corr*sd(x.perp)*X_aux + x.perp*sd(X_aux) * sqrt(1 - corr^2)
        X_vec <- rep(X, m)
        
      } else {
        
        ## X: continuous, time-varying; aux: continuous, time-invariant
        
        X_aux <- rnorm(n, mean = mean_X_aux, sd = sd_X_aux)
        X_aux_vec <- rep(X_aux, m)
        X_vec <- rep(NA, length(X_aux_vec))
        for (k in 1:m){
          x <- rnorm(n, mean = 0, sd = 1)
          x.perp <- residuals(lm(x ~ X_aux))
          X_vec[(1 + (k - 1)*n):(k*n)] <- corr*sd(x.perp)*X_aux + x.perp*sd(X_aux) * sqrt(1 - corr^2)
        }
        X <- X_vec[1:n]
        
      }
      
    } else {
      
      ## X: continuous, time-varying; aux: continuous, time-varying
      
      X_aux <- lapply(c(1:m), function(x){rnorm(n, mean = mean_X_aux, sd = sd_X_aux)})
      X_aux_vec <- unlist(X_aux)
      X_vec <- rep(NA, length(X_aux_vec))
      for (k in 1:m){
        x <- rnorm(n, mean = 0, sd = 1)
        X_aux_temp <- X_aux[[k]]
        x.perp <- residuals(lm(x ~ X_aux_temp))
        X_vec[(1 + (k - 1)*n):(k*n)] <- corr*sd(x.perp)*X_aux_temp + x.perp*sd(X_aux_temp) * sqrt(1 - corr^2)
      }
      
      X <- X_vec[1:n]
      
    }  
    
  } else if ((X.cont == TRUE) && (aux.cont == FALSE)){
    
    ## Auxiliary variables
    mean_X_aux <- 3
    sd_X_aux <- 0.1
    
    if (aux.timevar == FALSE){
      
      if (X.timevar == FALSE){
        
        ## X: continuous, time-invariant; aux: discrete, time-invariant 
        X_aux <- rnorm(n, mean = mean_X_aux, sd = sd_X_aux)
        X_aux[X_aux > mean_X_aux] <- 4
        X_aux[X_aux != 4] <- 2
        X_aux_vec <- rep(X_aux, m)
        x <- rnorm(n, mean = 0, sd = 1)
        x.perp <- residuals(lm(x ~ factor(X_aux)))
        X <- corr*sd(x.perp)*X_aux + x.perp*sd(X_aux) * sqrt(1 - corr^2)
        X_vec <- rep(X, m)
        
      } else {
        
        ## not yet modified
        X_aux <- rnorm(n, mean = mean_X_aux, sd = sd_X_aux)
        X_aux_vec <- rep(X_aux, m)
        X_vec <- rep(NA, length(X_aux_vec))
        for (k in 1:m){
          x <- rnorm(n, mean = 0, sd = 1)
          x.perp <- residuals(lm(x ~ X_aux))
          X_vec[(1 + (k - 1)*n):(k*n)] <- corr*sd(x.perp)*X_aux + x.perp*sd(X_aux) * sqrt(1 - corr^2)
        }
        X <- X_vec[1:n]
      }
      
    } else {
      
      ## X: continuous, time-varying; aux: discrete, time-varying 
      X_aux <- lapply(c(1:m), function(x){rnorm(n, mean = mean_X_aux, sd = sd_X_aux)})
      X_aux_vec <- unlist(X_aux)
      X_aux_vec[X_aux_vec > mean_X_aux] <- 4
      X_aux_vec[X_aux_vec != 4] <- 2
      X_vec <- rep(NA, length(X_aux_vec))
      for (k in 1:m){
        x <- rnorm(n, mean = 0, sd = 1)
        x.perp <- residuals(lm(x ~ factor(X_aux_vec[(1 + (k - 1)*n):(k*n)])))
        X_vec[(1 + (k - 1)*n):(k*n)] <- corr*sd(x.perp)*X_aux_vec[(1 + (k - 1)*n):(k*n)] + 
          x.perp*sd(X_aux_vec[(1 + (k - 1)*n):(k*n)]) * sqrt(1 - corr^2)
      }
      X <- X_vec[1:n]
    }
    
  } else if ((X.cont == FALSE) && (aux.cont == FALSE)){
    
    ## Auxiliary variables
    mean_X_aux <- 3
    sd_X_aux <- 0.1
    
    if (aux.timevar == FALSE){
      
      if (X.timevar == FALSE){
        
        X_aux <- rnorm(n, mean = mean_X_aux, sd = sd_X_aux)
        X_aux[X_aux > mean_X_aux] <- 4
        X_aux[X_aux != 4] <- 2
        X_aux_vec <- rep(X_aux, m)
        x <- rnorm(n, mean = 0, sd = 1)
        x.perp <- residuals(lm(x ~ factor(X_aux)))
        X <- corr*sd(x.perp)*X_aux + x.perp*sd(X_aux) * sqrt(1 - corr^2)
        thresX <- round(quantile(X)[4])
        X[X > thresX] <- 4
        X[X != 4] <- 2
        X_vec <- rep(X, m)
        
      } else {
        ## not yet modified
        X_aux <- rnorm(n, mean = mean_X_aux, sd = sd_X_aux)
        X_aux_vec <- rep(X_aux, m)
        X_vec <- rep(NA, length(X_aux_vec))
        for (k in 1:m){
          x <- rnorm(n, mean = 0, sd = 1)
          x.perp <- residuals(lm(x ~ X_aux))
          X_vec[(1 + (k - 1)*n):(k*n)] <- corr*sd(x.perp)*X_aux + x.perp*sd(X_aux) * sqrt(1 - corr^2)
        }
        X <- X_vec[1:n]
      }
      
    } else {
      ## not yet modified
      X_aux <- lapply(c(1:m), function(x){rnorm(n, mean = mean_X_aux, sd = sd_X_aux)})
      X_aux_vec <- unlist(X_aux)
      X_vec <- rep(NA, length(X_aux_vec))
      for (k in 1:m){
        x <- rnorm(n, mean = 0, sd = 1)
        X_aux_temp <- X_aux[[k]]
        x.perp <- residuals(lm(x ~ X_aux_temp))
        X_vec[(1 + (k - 1)*n):(k*n)] <- corr*sd(x.perp)*X_aux_temp + x.perp*sd(X_aux_temp) * sqrt(1 - corr^2)
      }
      
      X <- X_vec[1:n]
      
    }
    
  }
 
  return(list("X" = X,
              "X_vec" = X_vec,
              "X_aux" = X_aux,
              "X_aux_vec" = X_aux_vec))
   
}

