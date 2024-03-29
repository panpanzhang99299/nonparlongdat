##
## nonparlongdat: Nonparametric Methods for Longitudinal Data Analysis
## Copyright (C) 2024 Panpan Zhang and Sharon X. Xie
## Panpan Zhang <panpan.zhang@vumc.org>
##
## This file is part of the R package nonparlongdat.
##
## The R package wdnet is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package wdnet is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

#' @importFrom maxLik maxLik
#' @importFrom stats bw.SJ
NULL

## Nonparametric Estimation for Time-varying Covariates and Time-varying Auxiliary Variables

#' This function performs a nonparametric estimation with missing data in
#' time-varying covariates using time-varying auxilirary information.
#'
#' @param n.sub is a numeric input of the total number of subjects
#' @param m.visit is a numeric input of the total (maximum) number of visits for
#' each subject
#' @param data.y is a matrix/data frame of all longitudinal outcomes.
#' @param data.x is a matrix/data frame of time-varying covariates containing
#' missing values.
#' @param data.z is a matrix/data frame of fully observed covariates.
#' @param data.aux is a matrix/data frame of auxiliary variables.
#' @param para.ini is a vector of initial values of parameters.
#' @param cov.cont is a logistic argument given by \code{TRUE} if the missing
#' covariates are continuous; \code{FALSE}, otherwise.
#' @param aux.cont is a logistic argument given by \code{TRUE} if the auxiliary
#' variable is continuous; \code{FALSE}, otherwise.
#'
#' @return a list including parameter estimates, standard errors and the
#' convergence code
#'
#' @references To be added
#'
#' @note To be added
#'
#' @keywords internal
#'
#' @export


nonpar_X_timeinv_S_timeinv <- function(n.sub,
                                       m.visit,
                                       data.y,
                                       data.x,
                                       data.z,
                                       data.aux,
                                       para.ini,
                                       cov.cont = TRUE,
                                       aux.cont = TRUE) {
  n <- n.sub
  
  m <- m.visit
  
  Y_miss_wide <- matrix(data.y, ncol = m)
  
  X_miss_wide <- matrix(data.x, ncol = m)
  
  Z_wide <- matrix(data.z, ncol = m)
  
  X_aux_wide <- matrix(data.aux, ncol = m)
  
  
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
  
  if (aux.cont == TRUE) {
    ## Bandwidth
    H <- sapply(c(1:m), function(x) {
      stats::bw.SJ(X_aux_wide[V, x])
    })
    # H <- sapply(c(1:m), function(x){1.06*sd(X_aux_wide[V,x])*n^(-1/5)})
    
    loglik.nlme <- function(para) {
      theta1e <- para[1]
      
      theta2e <- para[2]
      
      theta3e <- para[3]
      
      lambdae <- para[4]
      
      sigmae <- para[5]
      
      ## The likelihood of the validation set
      theta.vec <- c(theta1e, theta2e, theta3e)
      
      Sigma.mat <- matrix(lambdae ^ 2, m, m) + diag(sigmae ^ 2, m)
      
      likcontri_valid <-
        loglikvalid (Y_miss_wide[V,],
                     X_miss_wide[V,],
                     Z_wide[V,],
                     theta.vec,
                     Sigma.mat)
      
      ## The likelihood of the nonvalidation set
      likcontri_nonv <-
        logliknonvalidXinvauxinv(
          Y_miss_wide[V,],
          Y_miss_wide[V_nonv,],
          X_miss_wide[V,],
          X_miss_wide[V_nonv,],
          Z_wide[V_nonv,],
          X_aux_wide[V,],
          X_aux_wide[V_nonv,],
          theta.vec,
          Sigma.mat,
          H,
          aux.cont
        )
      
      return(likcontri_valid + likcontri_nonv)
    }
    
  } else if (aux.cont == FALSE) {
    H <- rep(NA, m)
    
    loglik.nlme <- function(para) {
      theta1e <- para[1]
      
      theta2e <- para[2]
      
      theta3e <- para[3]
      
      lambdae <- para[4]
      
      sigmae <- para[5]
      
      ## The likelihood of the validation set
      theta.vec <- c(theta1e, theta2e, theta3e)
      
      Sigma.mat <- matrix(lambdae ^ 2, m, m) + diag(sigmae ^ 2, m)
      
      likcontri_valid <-
        loglikvalid (Y_miss_wide[V,],
                     X_miss_wide[V,],
                     Z_wide[V,],
                     theta.vec,
                     Sigma.mat)
      
      ## The likelihood of the nonvalidation set
      likcontri_nonv <-
        logliknonvalidXinvauxinv(
          Y_miss_wide[V,],
          Y_miss_wide[V_nonv,],
          X_miss_wide[V,],
          X_miss_wide[V_nonv,],
          Z_wide[V_nonv,],
          X_aux_wide[V,],
          X_aux_wide[V_nonv,],
          theta.vec,
          Sigma.mat,
          H,
          aux.cont
        )
      
      return(likcontri_valid + likcontri_nonv)
    }
    
  }
  
  para_mle <-
    maxLik::maxLik(loglik.nlme, start = para.ini, method = "BFGS")
  
  para_sum <- summary(para_mle)
  
  para_est <- para_mle$estimate
  
  para_se <- para_sum$estimate[, 2]
  
  para_code <- para_mle$code
  
  return(list(
    "est" = para_est,
    "se" = para_se,
    "converge" = para_code
  ))
  
}
