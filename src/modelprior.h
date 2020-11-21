#ifndef MODELPRIOR_H
#define MODELPRIOR_H

// Calculate area of hyper sphere
double area_s(double d){
  return(2*std::pow(arma::datum::pi,0.5*d+0.5)/tgamma(0.5*d+0.5));
}
// Simulate alpha for BPS_S kernel exiting hyper plane (use inversion)
double sim_alpha(double s){
  double sgn = 2.0*(R::runif(0,1) < 0.5) - 1;
  return(sqrt(1 - std::pow(R::runif(0,1),2/s))*sgn);
}
// Simulate alpha for BPS_S kernel exiting hyper plane (use Metropolis Hastings)
double q_density_S(double alpha, double s){
  return(0.5*std::abs(alpha)*s*std::pow(1-alpha*alpha, 0.5*(s-2)));
}
double MH_alpha_S(double s){
  double alpha_c = R::runif(-1,1);
  double q_c = q_density_S(alpha_c,s);
  double q_p, alpha_p;
  for( int k = 0; k < 15; k++ ){
    alpha_p = R::runif(-1,1);
    q_p = q_density_S(alpha_p, s);
    if( R::runif(0,1) < q_p/q_c){
      alpha_c = alpha_p;
      q_c = q_p;
    }
  }
  return(alpha_c);
}
// Simulate alpha for BPS_N kernel exiting hyper plane (use Metropolis Hastings)
double q_density(double alpha){
  return(2*std::abs(alpha)*exp(-0.5*alpha*alpha));
}
double MH_alpha_N(double s){
  double alpha_c = 2*(2.0*(R::runif(0,1) < 0.5) - 1);
  double q_c = q_density(alpha_c);
  double q_p, alpha_p;
  for( int k = 0; k < 10; k++ ){
    alpha_p = alpha_c + R::rnorm(0,1);
    q_p = q_density(alpha_p);
    if( R::runif(0,1) < q_p/q_c){
      alpha_c = alpha_p;
      q_c = q_p;
    }
  }
  return(alpha_c);
}

// Inital simple Idependent Gaussian prior
arma::vec zigzag_MP_GaussIID(arma::vec theta, double rjval, double ppi, double prior_sigma2) {
  arma::uvec inds_off = arma::find(arma::abs(theta) > 1e-10);
  arma::uvec inds_on = arma::find(arma::abs(theta) < 1e-10);
  int p = inds_on.size();
  arma::vec times(p);
  double alpha;

  for(int i =0; i < p; i++){
    alpha = rjval*ppi/(1-ppi)/sqrt(prior_sigma2*2*arma::datum::pi);
    times(i) = R::rexp(1)/alpha;
  }
  return(times);
}

arma::vec bps_MP_GaussIID_S(arma::vec theta, double rjval, double ppi, double prior_sigma2) {
  arma::uvec inds_off = arma::find(arma::abs(theta) > 1e-10);
  arma::uvec inds_on = arma::find(arma::abs(theta) < 1e-10);

  int num_off = inds_off.size();
  int num_on = theta.size() - num_off;
  double Jac = 1.0;

  if( num_off > 1e-10 ){
    Jac = area_s(theta.size() - num_on -1)/area_s(theta.size() - num_on);
    Jac = Jac*2/(theta.size() - num_on);
  }
  arma::vec times(num_on);

  double alpha;

  for(int i=0; i < num_on; i++){
    alpha = -0.5*log(2*arma::datum::pi*prior_sigma2);
    alpha += log(rjval) + log(ppi)-log(1-ppi) + log(Jac);
    times(i) = R::rexp(1)/exp(alpha);
  }
  return(times);
}
arma::vec bps_MP_GaussIID_N(arma::vec theta,double rjval, double ppi, double prior_sigma2) {
  arma::uvec inds_off = arma::find(arma::abs(theta) > 1e-10);
  arma::uvec inds_on = arma::find(arma::abs(theta) < 1e-10);

  int num_off = inds_off.size();
  int num_on = theta.size() - num_off;
  double Jac = 1.0;

  if( num_off > 1e-10 ){
    Jac = 2.0/sqrt(2.0*arma::datum::pi);
  }

  arma::vec times(num_on);
  double alpha;
  for(int i=0; i < num_on; i++){
    alpha = -0.5*log(2*arma::datum::pi*prior_sigma2);
    alpha += log(rjval) + log(ppi)-log(1-ppi) + log(Jac);
    times(i) = R::rexp(1)/exp(alpha);
  }
  return(times);
}
#endif
