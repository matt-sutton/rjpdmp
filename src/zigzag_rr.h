#ifndef ZIGZAG_RR_H
#define ZIGZAG_RR_H
#include "inversion.h"
#include "modelprior.h"
#include "likelyhood.h"

//' zigzag_rr
//'
//' Applies the Reversible Jump ZigZag Sampler to a Robust Regression problem with dirac spike and slab prior.
//' Included variables are given an independent Gaussian prior with variance \code{prior_sigma2} and mean zero, variables are given prior probabilities of inclusion \code{ppi}.
//'
//' @param maxTime Maximum runtime (in Seconds) of the algorithm; will terminate the code after a given computation time or nmax iterations of the algorithm is reached.
//' @param dataX Matrix of all covariates where the i-th row corresponds to all p covariates x_{i,1}, ..., x_{i,p} of the i-th observation.
//' @param datay Vector of n observations of a continuous response variable y.
//' @param prior_sigma2 Double for the prior variance for included variables.
//' @param x0 Initial position of the regression parameter
//' @param theta0 Initial velocity for the sampler (Default has 1s on all components). This should be chosen with unit velocities on each component (regardless of sign).
//' @param rj_val Reversible jump parameter for the PDMP method. This value is fixed over all models and is interpreted as the probability to jump to a reduced model when a parameter hits zero.
//' @param ppi Double for the prior probability of inclusion (ppi) for each parameter.
//' @param nmax Maximum number of iterations (simulated events) of the algorithm; will stop the algorithm when this number of iterations of the method have occured. Default value is 1e6, lower values should be chosen for memory constraints if less iterations are desired.
//' @param burn Optional number of iterations to use for burnin. These are not stored so can be useful in memory intensive problems.
//' @return Returns a list with the following objects:
//' @return \code{times}: Vector of event times where ZigZag process switchs velocity or jumps models.
//' @return \code{positions}: Matrix of positions at which event times occur, these are not samples from the PDMP.
//' @return \code{theta}: Matrix of new velocities at event times.
//' @examples
//'
//' generate.rr.data <- function(beta, n, Sig, noise, interc = TRUE) {
//' p <- length(beta)-(interc == TRUE)
//' dataX <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=Sig)
//' if(interc) {dataX <- cbind(1, dataX)}
//' dataY <- rep(0, n)
//' dataY <- dataX %*% as.vector(beta)+rnorm(n, sd = sqrt(noise))
//' return(list(dataX = dataX, dataY = dataY))
//' }
//' p <- 3;
//' n<- 120
//' beta <- c(0.5,0.5, rep(0,p-1))
//' set.seed(1)
//' data <- generate.rr.data(beta,n,diag(1,p+1), noise = 2, interc = FALSE)
//' dataX <- data$dataX; dataY <- data$dataY
//'
//' set.seed(1)
//' ppi_val <- 1/4
//' res <- zigzag_rr(maxTime = 1, dataX = dataX, datay = dataY,
//'                  prior_sigma2 = 1e2, x0 = rep(0,p+1), theta0 = rep(0,p+1),
//'                  rj_val = 0.6, ppi = ppi_val, nmax = 1e5)
//'\dontrun{
//' plot_pdmp(res, coords = 1:3, inds = 1:1e3)
//'}
//'
//' @export
// [[Rcpp::export]]
List zigzag_rr(double maxTime, const arma::mat& dataX, const arma::vec& datay,
               double prior_sigma2, arma::vec x0, arma::vec theta0, double rj_val = 0.5, double ppi=0.5,
               int nmax = 1e6, int burn = -1){
  int mini = 1, p = x0.size(), nEvent= 1;
  double eps = 1e-10, t = 0, grad_val, upper, val, tau_val;

  // Later functionality will take these as arguments
  double m = 0.5, sigma2_small = 1, sigma2_large = 1e2;

  arma::mat ab_vals(p,2), sk_points(p,nmax), sk_theta(p,nmax);
  arma::vec sk_times(nmax), theta = theta0, taus(p), x = x0;

  taus.zeros(); ab_vals.zeros(); sk_times.zeros();
  if( burn < 0){
    burn = 0;
    sk_points.col(0) = x0; sk_theta.col(0) = theta;
  }

  arma::vec grad_vals(p), a_vals(p), b_vals(p);

  arma::wall_clock timer;
  timer.tic();

  // Simulate initial times
  arma::uvec inds_off_hp = arma::find(arma::abs(theta) > eps);
  arma::uvec inds_on_hp = arma::find(arma::abs(theta) < eps);

  arma::vec resid = datay - dataX*x;
  arma::vec Xtheta = dataX*theta;
  arma::vec Xitheta(Xtheta.size());

  grad_vals.elem(inds_off_hp) =
    dataX.cols(inds_off_hp).t()*get_grad_resid_rr(resid, m, sigma2_large, sigma2_small) +
    x.elem(inds_off_hp)/prior_sigma2;

  a_vals.elem(inds_off_hp) = theta.elem(inds_off_hp)%grad_vals.elem(inds_off_hp);
  for(unsigned int i =0; i< inds_off_hp.size(); i++){
    b_vals(inds_off_hp(i)) = sum(arma::abs(theta(inds_off_hp(i))*dataX.col(inds_off_hp(i))%(Xtheta))) +
      theta(inds_off_hp(i))*theta(inds_off_hp(i))/prior_sigma2;
  }

  taus.elem(inds_off_hp) = get_time_off_hp(x.elem(inds_off_hp),theta.elem(inds_off_hp),a_vals.elem(inds_off_hp), b_vals.elem(inds_off_hp));
  taus.elem(inds_on_hp) = zigzag_MP_GaussIID(theta, rj_val, ppi, prior_sigma2);

  while( nEvent < nmax + burn ){
    mini = taus.index_min();
    tau_val = (taus(mini) - t);
    x += tau_val*theta;
    t = taus(mini);
    a_vals += b_vals*tau_val;
    resid -= Xtheta*tau_val;

    if( std::abs(x(mini)) < eps){
      x(mini) = 0;
      if( std::abs(theta(mini)) < eps ){
        // Exiting the HP
        theta[mini] = 2.0*(R::runif(0,1) < 0.5) - 1;
        Xtheta += dataX.col(mini)*theta(mini);
        if( nEvent >= burn){
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
      } else {
        // Enter the HP
        // w.p. rj enter the HP
        if(R::runif(0,1) < rj_val){
          Xtheta -= dataX.col(mini)*theta(mini);
          theta[mini] = 0;
          if( nEvent >= burn){
            sk_times(nEvent-burn) = t;
            sk_points.col(nEvent-burn) = x;
            sk_theta.col(nEvent-burn) = theta;
          }
          nEvent++;
        }
      }
      inds_off_hp = arma::find(arma::abs(theta) > eps);
      inds_on_hp = arma::find(arma::abs(theta) < eps);
      for(unsigned int i =0; i< inds_off_hp.size(); i++){
        b_vals(inds_off_hp(i)) = sum(arma::abs(theta(inds_off_hp(i))*dataX.col(inds_off_hp(i))%(Xtheta))) +
          theta(inds_off_hp(i))*theta(inds_off_hp(i))/prior_sigma2;
      }
      grad_vals.elem(inds_off_hp) =
        dataX.cols(inds_off_hp).t()*get_grad_resid_rr(resid, m, sigma2_large, sigma2_small) +
        x.elem(inds_off_hp)/prior_sigma2;
      a_vals.elem(inds_off_hp) = theta.elem(inds_off_hp)%grad_vals.elem(inds_off_hp);

    } else {

      upper = std::max(0.0, a_vals(mini));
      grad_val = arma::dot(dataX.col(mini),get_grad_resid_rr(resid, m, sigma2_large, sigma2_small)) +
        x(mini)/prior_sigma2;

      val = theta(mini)*grad_val;
      if( R::runif(0,1) < val/upper){
        Xtheta -= 2*dataX.col(mini)*theta(mini);
        theta[mini] = -theta[mini];
        if( nEvent >= burn){
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
        val = -val;
        for(unsigned int i =0; i< inds_off_hp.size(); i++){
          b_vals(inds_off_hp(i)) = sum(arma::abs(theta(inds_off_hp(i))*dataX.col(inds_off_hp(i))%(Xtheta))) +
            theta(inds_off_hp(i))*theta(inds_off_hp(i))/prior_sigma2;
        }
      }
      a_vals(mini) = val;
    }
    if(timer.toc() > maxTime){
      if(nEvent < burn){
        Rcout << "Sampler still in burnin phase - set a longer runtime" << std::endl;
      } else {
        sk_points.shed_cols(nEvent-burn, nmax-1);
        sk_theta.shed_cols(nEvent-burn, nmax-1);
        sk_times.shed_rows(nEvent-burn, nmax-1);
      }
      break;
    }

    taus.elem(inds_off_hp) = t + get_time_off_hp(x.elem(inds_off_hp),theta.elem(inds_off_hp),a_vals.elem(inds_off_hp), b_vals.elem(inds_off_hp));
    taus.elem(inds_on_hp) = t + zigzag_MP_GaussIID(theta, rj_val, ppi, prior_sigma2);
  }
  sk_times -= sk_times(0);

  List ret ;
  ret["times"] = sk_times ;
  ret["positions"] = sk_points ;
  ret["theta"] = sk_theta ;
  return(ret) ;
}

#endif
