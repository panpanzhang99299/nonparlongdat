#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;

//' @param X Vector, a vector of numeric values
//' @param mean Vector, the mean of multivariate normal
//' @param Sigma Matrix, the variance-covariance matrix of multivariate normal 
//' @return Returns the logarithm of the likelihood of X with multivariate normal with mean \code{mean} and variance \code{Sigma}
//'
//' @keywords internal
//'
 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double logdmvnorm (arma::rowvec const &X, // Input
                   arma::rowvec const &mean, // Mean vector
                   arma::mat const &Sigma // Variance-Covariance Matrix
){
  
  int p = X.n_cols;
  
  arma::mat const rootSigmainv = arma::inv(arma::trimatu(arma::chol(Sigma))); // The inverse of the Cholesky decomposion of Sigma
  
  double const Sigmaconst = arma::sum(log(rootSigmainv.diag())); // log of 1/|Simga|^(1/2)
  
  double const piconst = -(double)p/2 * std::log(2*M_PI); // log of 1/(2*pi)^(n/2)
  
  arma::rowvec z = (X - mean) * rootSigmainv;
  
  double out = Sigmaconst + piconst - 0.5 * arma::dot(z, z);
  
  return out;
                          
}

//' @param Xi Vector, a vector of numeric values
//' @param Xj Vector, a vector of numeric values
//' @return Returns a vector indentifying valid \code{Xi} for \code{Xj}
//' 
//' @keywords internal
//'

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::rowvec ideniforj(arma::rowvec const &Xi, 
                       arma::rowvec const &Xj)
{
  
  arma::uvec temp1 = arma::find_finite(Xi);
  
  arma::vec temp2 = Xj.t();
  
  arma::vec temp3 = temp2.elem(temp1);
  
  arma::rowvec out = temp3.t();
  
  return out;
  
}

//' @param Xi Vector, a vector of numeric values
//' @param Sigma Matrix, a variance-covariance matrix
//' @return Returns a matrix of correct dimension based on the index of \code{Xi}
//' 
//' @keywords internal
//'

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat idensigmaforj(arma::rowvec const &Xi, 
                        arma::mat const &Sigma)
{
  
  arma::uvec temp1 = arma::find_finite(Xi);
  
  arma::mat out = Sigma.submat(temp1, temp1);
  
  return out;
  
}

//' @param Y Matrix, a matrix of outcomes
//' @param X Matrix, a matrix with missing covariates
//' @param Z Matrix, a matrix with fully-observed covariates
//' @param theta Vector, a vector of unknown parameters
//' @param Sigma Matrix, a matrix of variance-covariance structure
//' @return Returns the log-likelihood of the validation set
//' 
//' @keywords internal
//'

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double loglikvalid (arma::mat const &Y,
                    arma::mat const &X,
                    arma::mat const &Z,
                    arma::rowvec const &theta,
                    arma::mat const &Sigma
){
  
  int n = Y.n_rows;
  
  arma::rowvec logliki(n);
  
  for (int i = 0; i < n; i++){
    
    arma::uvec indpres = arma::find_finite(Y.row(i));
    
    int numpres = indpres.n_elem;
    
    arma::rowvec temp1;
    
    temp1.ones(numpres);
    
    arma::mat temp2 = arma::join_vert(temp1, ideniforj(Y.row(i), X.row(i)));
    
    arma::mat temp3 = arma::join_vert(temp2, ideniforj(Y.row(i), Z.row(i)));
    
    arma::rowvec tempmean = theta * temp3;
    
    arma::mat tempsigma = idensigmaforj(Y.row(i), Sigma);
    
    arma::rowvec tempdensity = ideniforj(Y.row(i), Y.row(i));
    
    logliki(i) = logdmvnorm(tempdensity, tempmean, tempsigma);

  }
  
  double out = arma::sum(logliki);
  
  return out;
  
}

//' @param Yval Matrix, a matrix of outcomes from the validation set
//' @param Ynonval Matrix, a matrix of outcomes from the non-validation set
//' @param Xval Matrix, a matrix of missing covariates from the validation set
//' @param Xnonval Matrix, a matrix of missing covariates from the non-validation set
//' @param Z Matrix, a matrix with fully-observed covariates
//' @param auxval Matrix, a matrix of auxiliary variables from the validation set
//' @param auxnonval Matrix, a matrix of auxiliary variables from the non-validation set
//' @param theta Vector, a vector of unknown parameters
//' @param Sigma Matrix, a matrix of variance-covariance structure
//' @param H Vector, a vector of bandwidths
//' @param auxcont, a logistic argument determining whether the auxiliary variables are continuous
//' @return Returns the log-likelihood contributed by non-validation set with time-invariant X and time-invariant auxliary variables
//' 
//' @keywords internal
//'

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double logliknonvalidXinvauxinv (arma::mat const &Yval,
                                 arma::mat const &Ynonval,
                                 arma::mat const &Xval,
                                 arma::mat const &Xnonval,
                                 arma::mat const &Z,
                                 arma::mat const &auxval,
                                 arma::mat const &auxnonval,
                                 arma::rowvec const &theta,
                                 arma::mat const &Sigma,
                                 arma::rowvec const &H,
                                 bool &auxcont
){
  
  int nval = Yval.n_rows;
  
  int nonval = Ynonval.n_rows;
  
  int p = Yval.n_cols;
  
  arma::rowvec zeromean(p);
  
  zeromean.fill(0);
  
  arma::mat stdvarcov;
  
  stdvarcov.eye(p, p);
  
  arma::rowvec loglikj(nonval);
  
  for (int j = 0; j < nonval; j++){
    
    arma::rowvec kernelest(nval);
    
    for (int i = 0; i < nval; i++){
      
      if (auxcont == TRUE){
        
        double temp = (auxnonval(j, 0) - auxval(i, 0))/H(0);
        
        kernelest(i) = arma::normpdf(temp);
        
      } else if (auxcont == FALSE){
      
      if (auxnonval(j, 0) == auxval(i, 0)){
        
        kernelest(i) = 1;
        
      } else {
        
        kernelest(i) = 0; 
        
      }
      
      double matchingnum = arma::sum(kernelest);
      
      if (matchingnum < 10) {
        
        warning("Number of matching is less than 10!");
        
        }
      
      }
      
    }
    
    arma::uvec indpres = arma::find_finite(Ynonval.row(j));
    
    int numpres = indpres.n_elem;
    
    arma::rowvec YgivenXZ(nval);
    
    for (int i = 0; i < nval; i++){
      
      arma::rowvec temp1;
      
      temp1.ones(numpres);
      
      arma::mat temp2 = arma::join_vert(temp1, ideniforj(Ynonval.row(j), Xval.row(i)));
      
      arma::mat temp3 = arma::join_vert(temp2, ideniforj(Ynonval.row(j), Z.row(j)));
      
      arma::rowvec tempmean = theta * temp3;
      
      arma::mat tempsigma = idensigmaforj(Ynonval.row(j), Sigma);
      
      arma::rowvec tempdensity = ideniforj(Ynonval.row(j), Ynonval.row(j));
      
      YgivenXZ(i) = exp(logdmvnorm(tempdensity, tempmean, tempsigma));
      
    }
    
    loglikj(j) = log(arma::sum(arma::dot(YgivenXZ, kernelest))/arma::sum(kernelest));
      
    }
    
    double out = arma::sum(loglikj);
  
    return out;
  
}

//' @param Yval Matrix, a matrix of outcomes from the validation set
//' @param Ynonval Matrix, a matrix of outcomes from the non-validation set
//' @param Xval Matrix, a matrix of missing covariates from the validation set
//' @param Xnonval Matrix, a matrix of missing covariates from the non-validation set
//' @param Z Matrix, a matrix with fully-observed covariates
//' @param auxval Matrix, a matrix of auxiliary variables from the validation set
//' @param auxnonval Matrix, a matrix of auxiliary variables from the non-validation set
//' @param theta Vector, a vector of unknown parameters
//' @param Sigma Matrix, a matrix of variance-covariance structure
//' @param H Vector, a vector of bandwidths
//' @return Returns the log-likelihood contributed by non-validation set with time-varying X and time-invariant auxliary variables
//' 
//' @keywords internal
//'

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double logliknonvalidXvaryauxinv (arma::mat const &Yval,
                                  arma::mat const &Ynonval,
                                  arma::mat const &Xval,
                                  arma::mat const &Xnonval,
                                  arma::mat const &Z,
                                  arma::mat const &auxval,
                                  arma::mat const &auxnonval,
                                  arma::rowvec const &theta,
                                  arma::mat const &Sigma,
                                  arma::rowvec const &H
){
  
  int nval = Yval.n_rows;
  
  int nonval = Ynonval.n_rows;
  
  int p = Yval.n_cols;
  
  arma::rowvec zeromean(p);
  
  zeromean.fill(0);
  
  arma::mat stdvarcov;
  
  stdvarcov.eye(p, p);
  
  arma::rowvec loglikj(nonval);
  
  for (int j = 0; j < nonval; j++){
    
    arma::rowvec kernelinput(nval);
    
    for (int i = 0; i < nval; i++){
      
      kernelinput(i) = (auxnonval(j, 0) - auxval(i, 0))/H(0);
      
    }
    
    arma::rowvec kernelest = arma::normpdf(kernelinput);
    
    arma::uvec indpres = arma::find_finite(Ynonval.row(j));
    
    int numpres = indpres.n_elem;
    
    arma::rowvec YgivenXZ(nval);
    
    for (int i = 0; i < nval; i++){
      
      arma::rowvec temp1;
      
      temp1.ones(numpres);
      
      arma::rowvec temp2i = ideniforj(Ynonval.row(j), Xval.row(i));
      
      arma::rowvec temp2j = ideniforj(Ynonval.row(j), Xnonval.row(j));
      
      arma::uvec temp2missingind = arma::find_nonfinite(temp2j);
      
      temp2j(temp2missingind) = temp2i(temp2missingind);
        
      arma::mat temp2 = arma::join_vert(temp1, temp2j);
      
      arma::mat temp3 = arma::join_vert(temp2, ideniforj(Ynonval.row(j), Z.row(j)));
      
      arma::rowvec tempmean = theta * temp3;
      
      arma::mat tempsigma = idensigmaforj(Ynonval.row(j), Sigma);
      
      arma::rowvec tempdensity = ideniforj(Ynonval.row(j), Ynonval.row(j));
      
      YgivenXZ(i) = exp(logdmvnorm(tempdensity, tempmean, tempsigma));
      
      bool YgivenXZnan = std::isnan(YgivenXZ(i));
      
      if (YgivenXZnan == TRUE){
        
        YgivenXZ(i) = 0;
        
        kernelest(i) = 0;
        
      }
      
    }
    
    loglikj(j) = log(arma::sum(arma::dot(YgivenXZ, kernelest))/arma::sum(kernelest));
    
  }
  
  double out = arma::sum(loglikj);
  
  return out;
  
}

//' @param Yval Matrix, a matrix of outcomes from the validation set
//' @param Ynonval Matrix, a matrix of outcomes from the non-validation set
//' @param Xval Matrix, a matrix of missing covariates from the validation set
//' @param Xnonval Matrix, a matrix of missing covariates from the non-validation set
//' @param Z Matrix, a matrix with fully-observed covariates
//' @param auxval Matrix, a matrix of auxiliary variables from the validation set
//' @param auxnonval Matrix, a matrix of auxiliary variables from the non-validation set
//' @param theta Vector, a vector of unknown parameters
//' @param Sigma Matrix, a matrix of variance-covariance structure
//' @param H Vector, a vector of bandwidths
//' @param auxcont, a logistic argument determining whether the auxiliary variables are continuous
//' @return Returns the log-likelihood contributed by non-validation set with time-varying X and time-varying auxliary variables
//' 
//' @keywords internal
//'

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double logliknonvalidXvaryauxvary (arma::mat const &Yval,
                                   arma::mat const &Ynonval,
                                   arma::mat const &Xval,
                                   arma::mat const &Xnonval,
                                   arma::mat const &Z,
                                   arma::mat const &auxval,
                                   arma::mat const &auxnonval,
                                   arma::rowvec const &theta,
                                   arma::mat const &Sigma,
                                   arma::rowvec const &H,
                                   bool &auxcont
){
  
  int nval = Yval.n_rows;
  
  int nonval = Ynonval.n_rows;
  
  int p = Yval.n_cols;
  
  arma::rowvec zeromean(p);
  
  zeromean.fill(0);
  
  arma::mat stdvarcov;
  
  stdvarcov.eye(p, p);
  
  arma::rowvec loglikj(nonval);
  
  for (int j = 0; j < nonval; j++){
    
    arma::rowvec kernelest(nval);
    
    for (int i = 0; i < nval; i++){
      
      if (auxcont == TRUE){
        
        arma::rowvec tempauxj = ideniforj(Ynonval.row(j), auxnonval.row(j));
        
        arma::rowvec tempauxi = ideniforj(Ynonval.row(j), auxval.row(i));
        
        arma::rowvec tempauxH = ideniforj(Ynonval.row(j), H);
        
        arma::rowvec tempdensity = (tempauxj - tempauxi)/tempauxH;
        
        arma::rowvec tempmean = ideniforj(Ynonval.row(j), zeromean);
        
        arma::mat tempsigma = idensigmaforj(Ynonval.row(j), stdvarcov);
   
        kernelest(i) = exp(logdmvnorm(tempdensity, tempmean, tempsigma));
   
      } else if (auxcont == FALSE){
        
        bool tempbool = arma::approx_equal(auxnonval.row(j), auxval.row(i), "absdiff", 0.001);
        
        if (tempbool == TRUE){
          
          kernelest(i) = 1;
          
        } else {
          
          kernelest(i) = 0; 
          
        }
        
        double matchingnum = arma::sum(kernelest);
        
        if (matchingnum < 10) {
          
          warning("Number of matching is less than 10!");
          
        }
        
      }
      
    }
    
    arma::uvec indpres = arma::find_finite(Ynonval.row(j));
    
    int numpres = indpres.n_elem;
    
    arma::rowvec YgivenXZ(nval);
    
    for (int i = 0; i < nval; i++){
      
      arma::rowvec temp1;
      
      temp1.ones(numpres);
      
      arma::rowvec temp2i = ideniforj(Ynonval.row(j), Xval.row(i));
      
      arma::rowvec temp2j = ideniforj(Ynonval.row(j), Xnonval.row(j));
      
      arma::uvec temp2missingind = arma::find_nonfinite(temp2j);
      
      temp2j(temp2missingind) = temp2i(temp2missingind);
      
      arma::mat temp2 = arma::join_vert(temp1, temp2j);
      
      arma::mat temp3 = arma::join_vert(temp2, ideniforj(Ynonval.row(j), Z.row(j)));
      
      arma::rowvec tempmean = theta * temp3;
      
      arma::mat tempsigma = idensigmaforj(Ynonval.row(j), Sigma);
      
      arma::rowvec tempdensity = ideniforj(Ynonval.row(j), Ynonval.row(j));
      
      YgivenXZ(i) = exp(logdmvnorm(tempdensity, tempmean, tempsigma));
      
      bool YgivenXZnan = std::isnan(YgivenXZ(i));
      
      if (YgivenXZnan == TRUE){
        
        YgivenXZ(i) = 0;
        
        kernelest(i) = 0;
        
      }
      
    }
    
    loglikj(j) = log(arma::sum(arma::dot(YgivenXZ, kernelest))/arma::sum(kernelest));
    
  }
  
  double out = arma::sum(loglikj);
  
  return out;
  
}
  