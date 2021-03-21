#ifndef BPS_RR_H
#define BPS_RR_H
#include "inversion.h"
#include "likelyhood.h"
#include "modelprior.h"

//' bps_s_rr
//'
//' Applies the Reversible Jump BPS Sampler with Velocities drawn Uniformly on the p-Sphere to a Robust Regression problem with dirac spike and slab prior.
//' Included variables are given an independent Gaussian prior with variance \code{prior_sigma2} and mean zero, variables are given prior probabilities of inclusion \code{ppi}.
//'
//' @param maxTime Maximum runtime (in Seconds) of the algorithm; will terminate the code after a given computation time or nmax iterations of the algorithm is reached.
//' @param dataX Matrix of all covariates where the i-th row corresponds to all p covariates x_{i,1}, ..., x_{i,p} of the i-th observation.
//' @param datay Vector of n observations of a continuous response variable y.
//' @param prior_sigma2 Double for the prior variance for included variables.
//' @param x0 Initial position of the regression parameter
//' @param theta0 Initial velocity for the sampler.
//' @param ref Refreshment rate for BPS.
//' @param rj_val Reversible jump parameter for the PDMP method. This value is fixed over all models and is interpreted as the probability to jump to a reduced model when a parameter hits zero.
//' @param ppi Double for the prior probability of inclusion (ppi) for each parameter.
//' @param nmax Maximum number of iterations (simulated events) of the algorithm; will stop the algorithm when this number of iterations of the method have occured. Default value is 1e6, lower values should be chosen for memory constraints if less iterations are desired.
//' @param burn Optional number of iterations to use for burn-in. These are not stored so can be useful in memory intensive problems.
//' @return Returns a list with the following objects:
//' @return \code{times}: Vector of event times where ZigZag process switches velocity or jumps models.
//' @return \code{positions}: Matrix of positions at which event times occur, these are not samples from the PDMP.
//' @return \code{theta}: Matrix of new velocities at event times.
//' @examples
//'
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
//'\dontrun{
//' set.seed(1)
//' ppi_val <- 1/4
//' res <- bps_s_rr(maxTime = 1, dataX = dataX, datay = dataY,
//'                  prior_sigma2 = 1e2, x0 = rep(0,p+1), theta0 = rep(0,p+1),
//'                  rj_val = 0.6, ppi = ppi_val, nmax = 1e5)
//'
//' plot_pdmp(res, coords = 1:3, inds = 1:1e3)
//'}
//' @export
// [[Rcpp::export]]
List bps_s_rr(double maxTime, const arma::mat& dataX, const arma::vec& datay,
                 double prior_sigma2,arma::vec x0, arma::vec theta0, double ref = 0.1,
                 double rj_val = 0.5, double ppi=0.5, int nmax = 1e6, int burn = -1){
  int p = x0.size(), nEvent= 1;
  double eps = 1e-10, epsbig = 0.5, t = 0, upper, tau_val, alpha;

  arma::mat sk_points(p,nmax), sk_theta(p,nmax);
  arma::vec sk_times(nmax), theta = theta0, taus(p+2), x = x0, en_variables(p), ab_vals(2);

  // Later functionality will take these as arguments
  double m = 0.5, sigma2_small = 1, sigma2_large = 1e2;

  taus.zeros(); ab_vals.zeros(); sk_times.zeros(), en_variables.zeros();
  arma::wall_clock timer;
  timer.tic();

  // Simulate initial times
  // Assume linear bounds
  arma::uvec inds_off_hp = arma::find(arma::abs(theta) > eps);
  arma::uvec inds_on_hp = arma::find(arma::abs(theta) < eps);

  arma::vec resid = datay - dataX*x;
  arma::vec Xtheta = dataX*theta;

  en_variables.elem(inds_off_hp).ones();
  if(sum(en_variables) > epsbig ){
    theta = theta0/sqrt(arma::dot(theta0,theta0));
  }
  if( burn < 0){
    burn = 0;
    sk_points.col(0) = x0; sk_theta.col(0) = theta;
  }
  arma::vec grad_vals(p);

  grad_vals.elem(inds_off_hp) =
    dataX.cols(inds_off_hp).t()*get_grad_resid_rr(resid, m, sigma2_large, sigma2_small) +
    x.elem(inds_off_hp)/prior_sigma2;

  arma::mat Hess = dataX.t()*dataX;
  Hess.diag() += 1.0/prior_sigma2;

  ab_vals(0) = arma::dot(grad_vals, theta);
  ab_vals(1) = arma::dot(Hess*theta, theta);

  taus.elem(inds_off_hp) = get_hit_times(x.elem(inds_off_hp),theta.elem(inds_off_hp));
  taus.elem(inds_on_hp) = bps_MP_GaussIID_S(theta, rj_val, ppi, prior_sigma2);
  taus(p) = linear_inv_t(ab_vals(0),ab_vals(1),R::runif(0,1));
  taus(p+1) = R::rexp(1)/ref;

  // Rcout << " as ";
  while( nEvent < nmax + burn){
    int mini = taus.index_min();
    tau_val = (taus(mini) - t);
    x += tau_val*theta;
    t = taus(mini);
    ab_vals(0) += ab_vals(1)*tau_val;
    resid -= Xtheta*tau_val;

    if( mini == p + 1){
      // Rcout << " 1 ";
      // Refresh
      if(sum(en_variables) > epsbig ){
        theta.elem(inds_off_hp).randn();
        theta.elem(inds_off_hp) = theta.elem(inds_off_hp)/sqrt(arma::dot(theta.elem(inds_off_hp),theta.elem(inds_off_hp)));
      }
      Xtheta = dataX*theta;
      if( nEvent >= burn){
        sk_times(nEvent-burn) = t;
        sk_points.col(nEvent-burn) = x;
        sk_theta.col(nEvent-burn) = theta;
      }
      nEvent++;

      grad_vals.elem(inds_off_hp) =
        dataX.cols(inds_off_hp).t()*get_grad_resid_rr(resid, m, sigma2_large, sigma2_small) +
        x.elem(inds_off_hp)/prior_sigma2;
    }
    if( mini == p){
      // Rcout << "2";
      // Flip operator
      upper = std::max(0.0, ab_vals(0));

      grad_vals.elem(inds_off_hp) =
        dataX.cols(inds_off_hp).t()*get_grad_resid_rr(resid, m, sigma2_large, sigma2_small) +
        x.elem(inds_off_hp)/prior_sigma2;

      if( R::runif(0,1) < arma::dot(grad_vals,theta)/upper){

        theta.elem(inds_off_hp) -=
          2*arma::dot(theta.elem(inds_off_hp), grad_vals.elem(inds_off_hp))/
            arma::dot(grad_vals.elem(inds_off_hp),grad_vals.elem(inds_off_hp))*grad_vals.elem(inds_off_hp);

        Xtheta = dataX*theta;
        if( nEvent >= burn){
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
      }
    }
    if( mini < p) {
      // Rcout << " 3 ";
      x(mini) = 0;
      inds_off_hp = arma::find(en_variables > epsbig);
      if( en_variables(mini) > epsbig ){
        // w.p. rj enter HP
        if(R::runif(0,1) < rj_val){
          theta(mini) = 0;
          en_variables(mini) = 0;

          inds_off_hp = arma::find(en_variables > epsbig);
          inds_on_hp = arma::find(en_variables < epsbig);
          if(sum(en_variables) > epsbig ){
            theta.elem(inds_off_hp) = theta.elem(inds_off_hp)/sqrt(arma::dot(theta.elem(inds_off_hp),theta.elem(inds_off_hp)));
          }
          Xtheta = dataX*theta;

          if( nEvent >= burn){
            sk_times(nEvent-burn) = t;
            sk_points.col(nEvent-burn) = x;
            sk_theta.col(nEvent-burn) = theta;
          }
          nEvent++;
        }
      } else {
        // if mini is on HP then exit the HP
        if( inds_off_hp.size()  < epsbig ){
          alpha = 2.0*(R::runif(0,1) < 0.5) - 1;
        } else {
          alpha = sim_alpha(inds_off_hp.size());
        }
        theta(mini) = alpha;
        theta.elem(inds_off_hp) = theta.elem(inds_off_hp)/sqrt(arma::dot(theta.elem(inds_off_hp),theta.elem(inds_off_hp)));

        Xtheta = dataX*theta;
        en_variables(mini) = 1;
        inds_off_hp = arma::find(en_variables > epsbig);
        inds_on_hp = arma::find(en_variables < epsbig);
        if( nEvent >= burn){
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
      }

      grad_vals.zeros();
      grad_vals.elem(inds_off_hp) =
        dataX.cols(inds_off_hp).t()*get_grad_resid_rr(resid, m, sigma2_large, sigma2_small) +
        x.elem(inds_off_hp)/prior_sigma2;
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

    ab_vals(0) = arma::dot(grad_vals.elem(inds_off_hp), theta.elem(inds_off_hp));
    ab_vals(1) = arma::dot(Hess(inds_off_hp,inds_off_hp)*theta.elem(inds_off_hp), theta.elem(inds_off_hp));

    taus.elem(inds_off_hp) = t + get_hit_times(x.elem(inds_off_hp),theta.elem(inds_off_hp));
    taus.elem(inds_on_hp) = t + bps_MP_GaussIID_S(theta, rj_val, ppi, prior_sigma2);
    taus(p) = t + linear_inv_t(ab_vals(0),ab_vals(1),R::runif(0,1));
    taus(p+1) = t + R::rexp(1)/ref;
  }
  sk_times -= sk_times(0);

  List ret ;
  ret["times"] = sk_times ;
  ret["positions"] = sk_points ;
  ret["theta"] = sk_theta ;
  return(ret) ;
}

//' bps_n_rr
//'
//' Applies the Reversible Jump BPS Sampler with Velocities drawn from the Normal distribution to a Robust Regression problem with dirac spike and slab prior.
//' Included variables are given an independent Gaussian prior with variance \code{prior_sigma2} and mean zero, variables are given prior probabilities of inclusion \code{ppi}.
//'
//' @param maxTime Maximum runtime (in Seconds) of the algorithm; will terminate the code after a given computation time or nmax iterations of the algorithm is reached.
//' @param dataX Matrix of all covariates where the i-th row corresponds to all p covariates x_{i,1}, ..., x_{i,p} of the i-th observation.
//' @param datay Vector of n observations of a continuous response variable y.
//' @param prior_sigma2 Double for the prior variance for included variables.
//' @param x0 Initial position of the regression parameter
//' @param theta0 Initial velocity for the sampler.
//' @param ref Refreshment rate for BPS.
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
//'\dontrun{
//' set.seed(1)
//' ppi_val <- 1/4
//' res <- bps_n_rr(maxTime = 1, dataX = dataX, datay = dataY,
//'                  prior_sigma2 = 1e2, x0 = rep(0,p+1), theta0 = rep(0,p+1),
//'                  rj_val = 0.6, ppi = ppi_val, nmax = 1e5, ref = 0.1, burn = -1)
//'
//' plot_pdmp(res, coords = 1:3, inds = 1:1e3)
//'}
//' @export
// [[Rcpp::export]]
List bps_n_rr(double maxTime, const arma::mat& dataX, const arma::vec& datay,
                 double prior_sigma2, arma::vec x0, arma::vec theta0, double ref = 0.1, double rj_val = 0.5,
                 double ppi=0.5, int nmax = 1e6, int burn = -1){

  int p = x0.size(), nEvent= 1;
  double eps = 1e-10, epsbig = 0.5, t = 0, upper, tau_val, alpha;

  // Later functionality will take these as arguments
  double m = 0.5, sigma2_small = 1, sigma2_large = 1e2;

  arma::mat sk_points(p,nmax), sk_theta(p,nmax);
  arma::vec sk_times(nmax), theta = theta0, taus(p+2), x = x0, en_variables(p), ab_vals(2);

  taus.zeros(); ab_vals.zeros(); sk_times.zeros(), en_variables.zeros();
  if( burn < 0){
    burn = 0;
    sk_points.col(0) = x0; sk_theta.col(0) = theta;
  }
  arma::wall_clock timer;
  timer.tic();

  // Simulate initial times
  // Assume linear bounds
  arma::uvec inds_off_hp = arma::find(arma::abs(theta) > eps);
  arma::uvec inds_on_hp = arma::find(arma::abs(theta) < eps);

  arma::vec resid = datay - dataX*x;
  arma::vec Xtheta = dataX*theta;

  en_variables.elem(inds_off_hp).ones();
  arma::vec grad_vals(p);

  grad_vals.elem(inds_off_hp) =
    dataX.cols(inds_off_hp).t()*get_grad_resid_rr(resid, m, sigma2_large, sigma2_small) +
    x.elem(inds_off_hp)/prior_sigma2;

  arma::mat Hess = dataX.t()*dataX;
  Hess.diag() += 1.0/prior_sigma2;

  ab_vals(0) = arma::dot(grad_vals, theta);
  ab_vals(1) = arma::dot(Hess*theta, theta);

  taus.elem(inds_off_hp) = get_hit_times(x.elem(inds_off_hp),theta.elem(inds_off_hp));
  taus.elem(inds_on_hp) = bps_MP_GaussIID_N(theta, rj_val, ppi, prior_sigma2);
  taus(p) = linear_inv_t(ab_vals(0),ab_vals(1),R::runif(0,1));
  taus(p+1) = R::rexp(1)/ref;

  // Rcout << " as ";
  while( nEvent < nmax + burn){
    int mini = taus.index_min();
    tau_val = (taus(mini) - t);
    x += tau_val*theta;
    t = taus(mini);
    ab_vals(0) += ab_vals(1)*tau_val;
    resid -= Xtheta*tau_val;

    if( mini == p + 1){
      // Refresh
      theta.elem(arma::find(en_variables > epsbig)).randn();
      Xtheta = dataX*theta;
      if( nEvent >= burn){
        sk_times(nEvent-burn) = t;
        sk_points.col(nEvent-burn) = x;
        sk_theta.col(nEvent-burn) = theta;
      }
      nEvent++;

      grad_vals.elem(inds_off_hp) =
        dataX.cols(inds_off_hp).t()*get_grad_resid_rr(resid, m, sigma2_large, sigma2_small) +
        x.elem(inds_off_hp)/prior_sigma2;
    }
    if( mini == p){
      // Rcout << "2";
      // Flip operator
      upper = std::max(0.0, ab_vals(0));

      grad_vals.elem(inds_off_hp) =
        dataX.cols(inds_off_hp).t()*get_grad_resid_rr(resid, m, sigma2_large, sigma2_small) +
        x.elem(inds_off_hp)/prior_sigma2;

      if( R::runif(0,1) < arma::dot(grad_vals,theta)/upper){

        theta.elem(inds_off_hp) -=
          2*arma::dot(theta.elem(inds_off_hp), grad_vals.elem(inds_off_hp))/
            arma::dot(grad_vals.elem(inds_off_hp),grad_vals.elem(inds_off_hp))*grad_vals.elem(inds_off_hp);

        Xtheta = dataX*theta;
        if( nEvent >= burn){
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
      }
    }
    if( mini < p) {
      // Rcout << "3";
      x(mini) = 0;
      inds_off_hp = arma::find(en_variables > epsbig);
      if( en_variables(mini) > epsbig ){
        // w.p. rj enter HP
        if(R::runif(0,1) < rj_val){
          theta[mini] = 0;
          Xtheta = dataX*theta;

          en_variables(mini) = 0;
          if( nEvent >= burn){
            sk_times(nEvent-burn) = t;
            sk_points.col(nEvent-burn) = x;
            sk_theta.col(nEvent-burn) = theta;
          }
          nEvent++;
        }
      } else {
        // if mini is on HP then exit the HP
        if( inds_off_hp.size()  < epsbig ){
          alpha = 2.0*(R::runif(0,1) < 0.5) - 1;
        } else {
          alpha = MH_alpha_N(inds_off_hp.size());
        }
        theta[mini] = alpha;

        Xtheta = dataX*theta;
        en_variables(mini) = 1;
        if( nEvent >= burn){
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
      }

      inds_off_hp = arma::find(en_variables > epsbig);
      inds_on_hp = arma::find(en_variables < epsbig);

      grad_vals.zeros();
      grad_vals.elem(inds_off_hp) =
        dataX.cols(inds_off_hp).t()*get_grad_resid_rr(resid, m, sigma2_large, sigma2_small) +
        x.elem(inds_off_hp)/prior_sigma2;
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

    ab_vals(0) = arma::dot(grad_vals.elem(inds_off_hp), theta.elem(inds_off_hp));
    ab_vals(1) = arma::dot(Hess(inds_off_hp,inds_off_hp)*theta.elem(inds_off_hp), theta.elem(inds_off_hp));

    taus.elem(inds_off_hp) = t + get_hit_times(x.elem(inds_off_hp),theta.elem(inds_off_hp));
    taus.elem(inds_on_hp) = t + bps_MP_GaussIID_N(theta, rj_val, ppi, prior_sigma2);
    taus(p) = t + linear_inv_t(ab_vals(0),ab_vals(1),R::runif(0,1));
    taus(p+1) = t + R::rexp(1)/ref;
  }
  sk_times -= sk_times(0);

  List ret ;
  ret["times"] = sk_times ;
  ret["positions"] = sk_points ;
  ret["theta"] = sk_theta ;
  return(ret) ;
}
#endif
