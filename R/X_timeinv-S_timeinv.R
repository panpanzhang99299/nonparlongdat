nonpara_X_timeinv_S_timeinv <- function(lme_df_miss,
                                        para_ini,
                                        X.cont = TRUE, 
                                        aux.cont = TRUE){
  
  library(Rcpp)
  
  sourceCpp("kernelest.cpp")
  
  Y_miss_wide <- matrix(lme_df_miss$resp, ncol = m)
  
  X_miss_wide <- matrix(lme_df_miss$covX, ncol = m)

  Z_wide <- matrix(lme_df_miss$covZ, ncol = m)
  
  X_aux_wide <- matrix(lme_df_miss$aux, ncol = m)
  
  ## Missing value indicator
  temp <- Y_miss_wide + X_miss_wide
  temp[is.na(temp)] = 0
  
  ## Identify the validation set
  V <- which(rowSums(temp) > 0)
  
  ## Number of the subjects in the validation set
  v <- length(V)
  
  ## Identify the completely (or partially) missing set
  V_nonv <- setdiff(c(1:n), V)
  
  ## Number of the subjects in the non-validation set
  v_nonv <- length(V_nonv)
  
  ##------------------------------------------------------------------------------------------------------------
  ########################################
  #### Compute the full log liklihood ####
  ########################################
  
  if (aux.cont == TRUE){
    
    ## Bandwidth
    H <- sapply(c(1:m), function(x){bw.SJ(X_aux_wide[V, x])})
    ## H <- sapply(c(1:m), function(x){1.06*sd(X_aux_wide[V,x])*n^(-1/5)})
    
    loglik.nlme <- function(para){
      
      theta1e <- para[1]
      
      theta2e <- para[2]
      
      theta3e <- para[3]
      
      lambdae <- lambda
      
      sigmae <- sigma
      
      ## The likelihood of the validation set
      
      theta.vec <- c(theta1e, theta2e, theta3e)
      
      Sigma.mat <- matrix(lambdae^2, m, m) + diag(sigmae^2, m)
      
      likcontri_valid <- loglikvalid (Y_miss_wide[V,], X_miss_wide[V,], Z_wide[V,], theta.vec, Sigma.mat)
      
      ## The likelihood of the nonvalidation set
      
      likcontri_nonv <- logliknonvalidXinvauxinv(Y_miss_wide[V,], Y_miss_wide[V_nonv,], X_miss_wide[V,], X_miss_wide[V_nonv,], Z_wide[V_nonv,],
                                                 X_aux_wide[V,], X_aux_wide[V_nonv,], theta.vec, Sigma.mat, H, aux.cont)
      
      return(likcontri_valid + likcontri_nonv)
    }
    
  } else if (aux.cont == FALSE){
    
    H <- rep(NA, m)
    
    loglik.nlme <- function(para){
      
      theta1e <- para[1]
      
      theta2e <- para[2]
      
      theta3e <- para[3]
      
      lambdae <- lambda
      
      sigmae <- sigma
      
      ## The likelihood of the validation set
      
      theta.vec <- c(theta1e, theta2e, theta3e)
      
      Sigma.mat <- matrix(lambdae^2, m, m) + diag(sigmae^2, m)
      
      likcontri_valid <- loglikvalid (Y_miss_wide[V,], X_miss_wide[V,], Z_wide[V,], theta.vec, Sigma.mat)
      
      ## The likelihood of the nonvalidation set
      
      likcontri_nonv <- logliknonvalidXinvauxinv(Y_miss_wide[V,], Y_miss_wide[V_nonv,], X_miss_wide[V,], X_miss_wide[V_nonv,], Z_wide[V_nonv,],
                                                 X_aux_wide[V,], X_aux_wide[V_nonv,], theta.vec, Sigma.mat, H, aux.cont)
      
      return(likcontri_valid + likcontri_nonv)
    }
    
  }  
  
  para_mle <- maxLik(loglik.nlme, start = para_ini, method = "BFGS")
  
  para_sum <- summary(para_mle)
  
  para_est<- para_mle$estimate
  
  para_se <- para_sum$estimate[,2]
  
  para_code <- para_mle$code
  
  return(list("est" = para_est, "se" = para_se, "converge" = para_code))
  
}

