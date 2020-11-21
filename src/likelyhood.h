#ifndef LIKELYHOOD_H
#define LIKELYHOOD_H

// LOGISTIC //
// Log likelihood (useful for finding control variate) //
double log_logit(arma::mat X, arma::vec y, arma::vec x){
  arma::vec z = X*x;
  return(  sum( log(1+exp(z)) - y%z)) ;
}

// Gradient of the Logistic Regression likelihood //
arma::vec get_grad_logit(const arma::mat& X, arma::vec z,const arma::vec& y, arma::uvec partials){
  int grad_size = partials.size();
  arma::vec grad(grad_size);
  for( int i = 0; i < grad_size; i++){
    grad(i) = sum(X.col(partials(i))/(1+exp(-z)) - y%X.col(partials(i)));
  }
  return(grad);
}

// Subsampled Gradient of the Logistic Regression likelihood //
arma::vec get_grad_logit_ss(const arma::mat& X, arma::vec x, const arma::vec& y,
                                  arma::uvec partials, arma::vec cvref, bool use_cv, arma::vec gradcv){
  int grad_size = partials.size(), n_obs = y.size();
  arma::vec grad(grad_size);
  unsigned int J = floor(n_obs*R::runif(0,1));
  double z = arma::dot(X.row(J),x);
  if(use_cv){
    double z_cv = arma::dot(X.row(J),cvref);
    for( int i = 0; i < grad_size; i++){
      grad(i) = gradcv(partials(i)) + n_obs*X(J,partials(i))*(1.0/(1.0+exp(-z)) -
        1.0/(1.0+exp(-z_cv)));
    }
  } else {
    for( int i = 0; i < grad_size; i++){
      grad(i) = n_obs*(X(J,partials(i))/(1+exp(-z)) - y(J)*X(J,partials(i)));
    }
  }
  return(grad);
}


// ROBUST REGRESSION //

// Aux function
arma::vec log_add_exp(arma::vec a, arma::vec b){
  arma::vec c = arma::max(a, b);
  return(c + log(exp(a-c) + exp(b-c)));
}

// Gradient of residual from likelihood g'(e(t))//
arma::vec get_grad_resid_rr(arma::vec resid, double m, double sig12, double sig22){

  double sig16 = std::pow(sig12, 3); double sig26 = std::pow(sig22, 3);
  arma::vec resid2 = resid%resid;

  // resid_j ~ (1-m)N(0,sig1) + mN(0,sig2)
  arma::vec temp1 = log_add_exp(log(1-m) - 0.5*log(sig12*2*arma::datum::pi) - resid2/(2*sig12),
                                log(m) - 0.5*log(sig22*2*arma::datum::pi) - resid2/(2*sig22));

  arma::vec temp2 = log_add_exp(log(1-m) - 0.5*log(sig16*2*arma::datum::pi) - resid2/(2*sig12),
                                log(m) - 0.5*log(sig26*2*arma::datum::pi) - resid2/(2*sig22));
  arma::vec temp3 = resid%exp(temp2 - temp1);
  return(-1.0*temp3);
}

#endif
