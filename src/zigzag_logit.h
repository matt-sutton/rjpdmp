#ifndef ZIGZAG_LOGIT_H
#define ZIGZAG_LOGIT_H

#include "inversion.h"
#include "modelprior.h"
#include "likelyhood.h"


//' zigzag_logit
//'
//' Applies the Reversible Jump ZigZag Sampler to a Logistic regression with dirac spike and slab distribution, as detailed in Reversible Jump PDMP Samplers for Variable Selection, 2020.
//' For included variables an independent Gaussian prior is assumed with variance \code{prior_sigma2} and mean zero, variables are given prior probabilities of inclusion \code{ppi}.
//'
//' @param maxTime Maximum runtime (in Seconds) of the algorithm; will terminate the code after a given computation time or nmax iterations of the algorithm is reached.
//' @param dataX Matrix of all covariates where the i-th row corresponds to all p covariates x_{i,1}, ..., x_{i,p} of the i-th observation.
//' @param datay Vector of n observations of a {0, 1}-valued variable y.
//' @param prior_sigma2 Double for the prior variance for included variables.
//' @param x0 Initial position of the regression parameter
//' @param theta0 Initial velocity for the sampler (Default has 1s on all components). This should be chosen with unit velocities on each component (regardless of sign).
//' @param rj_val Reversible jump parameter for the PDMP method. This value is fixed over all models and is interpreted as the probability to jump to a reduced model when a parameter hits zero.
//' @param ppi Double for the prior probability of inclusion (ppi) for each parameter.
//' @param nmax Maximum number of iterations (simulated events) of the algorithm; will stop the algorithm when this number of iterations of the method have occured. Default value is 10^6, lower values should be chosen for memory constraints if less iterations are desired.
//' @param burn Optional number of iterations to use for burnin. These are not stored so can be useful in memory intensive problems.
//' @return Returns a list with the following objects:
//' @return \code{times}: Vector of event times where ZigZag process switchs velocity or jumps models.
//' @return \code{positions}: Matrix of positions at which event times occur, these are not samples from the PDMP.
//' @return \code{theta}: Matrix of new velocities at event times.
//' @examples
//'
//' generate.logistic.data <- function(beta, n.obs, Sig) {
//' p <- length(beta)
//' dataX <- MASS::mvrnorm(n=n.obs,mu=rep(0,p),Sigma=Sig)
//' vals <- dataX %*% as.vector(beta)
//' generateY <- function(p) { rbinom(1, 1, p)}
//' dataY <- sapply(1/(1 + exp(-vals)), generateY)
//' return(list(dataX = dataX, dataY = dataY))
//' }
//'
//' n <- 15
//' p <- 25
//' beta <- c(1, rep(0, p-1))
//' Siginv <- diag(1,p,p)
//' Siginv[1,2] <- Siginv[2,1] <- 0.9
//' set.seed(1)
//' data <- generate.logistic.data(beta, n, solve(Siginv))
//' ppi <- 2/p
//'
//' zigzag_fit <- zigzag_logit(maxTime = 5, dataX = data$dataX, datay = data$dataY,
//'                            prior_sigma2 = 10,theta0 = rep(0, p), x0 = rep(0, p), rj_val = 0.6,
//'                            ppi = ppi)
//'
//' gibbs_fit <- gibbs_logit(maxTime = 5, dataX = data$dataX, datay = data$dataY,
//'                          prior_sigma2 = 10,beta = rep(0,p), gamma = rep(0,p),
//'                          ppi = ppi)
//'
//' plot_pdmp(zigzag_fit, coords = 1:2, inds = 1:10^4,burn = .1, nsamples = 5*1e4,
//'            mcmc_samples = t(gibbs_fit$beta*gibbs_fit$gamma))
//'
//' @export
// [[Rcpp::export]]
List zigzag_logit(double maxTime, const arma::mat& dataX, const arma::vec& datay,
                    double prior_sigma2, arma::vec x0, arma::vec theta0,
                    double rj_val = 0.6, double ppi=0.5, int nmax = 10^6, int burn = -1){
  int mini = 1, p = x0.size(), nEvent= 1;
  double eps = 1e-10, t = 0, upper, val, tau_val;

  arma::mat ab_vals(p,2), sk_points(p,nmax), sk_theta(p,nmax);
  arma::vec sk_times(nmax), theta = theta0, taus(p), x = x0, grad_val(1);

  taus.zeros(); ab_vals.zeros(); sk_times.zeros();
  if( burn < 0){
    burn = 0;
    sk_points.col(0) = x0; sk_theta.col(0) = theta;
  }

  arma::vec grad_vals(p), a_vals(p), b_vals(p);
  b_vals.zeros(); a_vals.zeros();

  arma::wall_clock timer;
  timer.tic();

  // Simulate initial times
  arma::uvec inds_off_hp = arma::find(arma::abs(theta) > eps);
  arma::uvec inds_on_hp = arma::find(arma::abs(theta) < eps);
  arma::uvec miniu(1);

  arma::vec Xtheta = dataX*theta;
  arma::vec z = dataX*x;

  grad_vals.elem(inds_off_hp) = get_grad_logit(dataX, z, datay, inds_off_hp) +
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
    z += Xtheta*tau_val;

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
      grad_vals.elem(inds_off_hp) = get_grad_logit(dataX, z, datay, inds_off_hp) +
        x.elem(inds_off_hp)/prior_sigma2;
      a_vals.elem(inds_off_hp) = theta.elem(inds_off_hp)%grad_vals.elem(inds_off_hp);

    } else {

      upper = std::max(0.0, a_vals(mini));
      miniu(0) = mini;
      grad_val = get_grad_logit(dataX, z, datay, miniu) + x(mini)/prior_sigma2;

      val = theta(mini)*grad_val(0);
      if( R::runif(0,1) < val/upper){
        Xtheta -= 2*dataX.col(mini)*theta(mini);
        theta(mini) = -theta(mini);
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

    taus.elem(inds_off_hp) = t + get_time_off_hp(x.elem(inds_off_hp),theta.elem(inds_off_hp),a_vals.elem(inds_off_hp),b_vals.elem(inds_off_hp));
    taus.elem(inds_on_hp) = t + zigzag_MP_GaussIID(theta, rj_val, ppi, prior_sigma2);
  }
  sk_times -= sk_times(0);

  List ret ;
  ret["times"] = sk_times ;
  ret["positions"] = sk_points ;
  ret["theta"] = sk_theta ;
  return(ret) ;
}

// Aux function to check if the control variate should be used //
bool check_cv(arma::uvec inds_off_hp_cv, arma::uvec inds_off_hp){
  if( arma::approx_equal(inds_off_hp_cv, inds_off_hp, "absdiff", 0.002) ){
    return(true);
  } else {
    return(false);
  }
}

//' zigzag_logit_ss
//'
//' Applies the Reversible Jump ZigZag Sampler with subsampling to a Logistic regression with dirac spike and slab distribution, as detailed in Reversible Jump PDMP Samplers for Variable Selection, 2020.
//' For included variables an independent Gaussian prior is assumed with variance \code{prior_sigma2} and mean zero, variables are given prior probabilities of inclusion \code{ppi}.
//'
//' @param maxTime Maximum runtime (in Seconds) of the algorithm; will terminate the code after a given computation time or nmax iterations of the algorithm is reached.
//' @param dataX Matrix of all covariates where the i-th row corresponds to all p covariates x_{i,1}, ..., x_{i,p} of the i-th observation.
//' @param datay Vector of n observations of a {0, 1}-valued variable y.
//' @param prior_sigma2 Double for the prior variance for included variables.
//' @param x0 Initial position of the regression parameter
//' @param theta0 Initial velocity for the sampler (Default has 1s on all components). This should be chosen with unit velocities on each component (regardless of sign).
//' @param cvref Control variate vector of dimension p for subsampling. If no control variate set to a vector of zeros.
//' @param rj_val Reversible jump parameter for the PDMP method. This value is fixed over all models and is interpreted as the probability to jump to a reduced model when a parameter hits zero.
//' @param ppi Double for the prior probability of inclusion (ppi) for each parameter.
//' @param nmax Maximum number of iterations (simulated events) of the algorithm; will stop the algorithm when this number of iterations of the method have occured. Default value is 10^6, lower values should be chosen for memory constraints if less iterations are desired.
//' @param burn Optional number of iterations to use for burnin. These are not stored so can be useful in memory intensive problems.
//' @return Returns a list with the following objects:
//' @return \code{times}: Vector of event times where ZigZag process switchs velocity or jumps models.
//' @return \code{positions}: Matrix of positions at which event times occur, these are not samples from the PDMP.
//' @return \code{theta}: Matrix of new velocities at event times.
//' @examples
//'
//' generate.logistic.data <- function(beta, n.obs, Sig) {
//' p <- length(beta)
//' dataX <- MASS::mvrnorm(n=n.obs,mu=rep(0,p),Sigma=Sig)
//' vals <- dataX %*% as.vector(beta)
//' generateY <- function(p) { rbinom(1, 1, p)}
//' dataY <- sapply(1/(1 + exp(-vals)), generateY)
//' return(list(dataX = dataX, dataY = dataY))
//' }
//'
//' n <- 15
//' p <- 25
//' beta <- c(1, rep(0, p-1))
//' Siginv <- diag(1,p,p)
//' Siginv[1,2] <- Siginv[2,1] <- 0.9
//' set.seed(1)
//' data <- generate.logistic.data(beta, n, solve(Siginv))
//' ppi <- 2/p
//'
//' zigzag_fit <- zigzag_logit(maxTime = 5, dataX = data$dataX,
//'                            datay = data$dataY, prior_sigma2 = 10,
//'                            theta0 = rep(0, p), x0 = rep(0, p), rj_val = 0.6,
//'                            ppi = ppi)
//'
//' zigzag_fit_s <- zigzag_logit_ss(maxTime = 5, dataX = data$dataX,
//'                                 datay = data$dataY,prior_sigma2 = 10,
//'                                 theta0 = rep(0, p), x0 = rep(0, p),
//'                                 rj_val = 0.6, cvref = c(1,rep(0,p-1)),
//'                                 ppi = ppi)
//'
//' gibbs_fit <- gibbs_logit(maxTime = 5, dataX = data$dataX, datay =data$dataY,
//'                          prior_sigma2 = 10,beta = rep(0,p), gamma =rep(0,p),
//'                          ppi = ppi)
//'
//' plot_pdmp_multiple(list(zigzag_fit,zigzag_fit_s), coords = 1:2, burn = .1,
//'                    inds = 1:length(zigzag_fit_s$times), nsamples = 5*1e4,
//'                    mcmc_samples = t(gibbs_fit$beta*gibbs_fit$gamma))
//'
//' @export
// [[Rcpp::export]]
List zigzag_logit_ss(double maxTime, const arma::mat& dataX, const arma::vec& datay,
                    double prior_sigma2, arma::vec x0, arma::vec theta0, arma::vec cvref,
                    double rj_val = 0.6, double ppi=0.5, int nmax = 10^6, int burn = -1){
  int mini = 1, p = x0.size(), nEvent= 1, n_obs = dataX.n_rows;
  double eps = 1e-10, t = 0, upper, val, tau_val;
  bool use_cv = false;

  arma::mat ab_vals(p,2), sk_points(p,nmax), sk_theta(p,nmax);
  arma::vec sk_times(nmax), theta = theta0, taus(p), x = x0, grad_val(1);
  taus.zeros(); ab_vals.zeros(); sk_times.zeros();
  if( burn < 0){
    burn = 0;
    sk_points.col(0) = x0; sk_theta.col(0) = theta;
  }

  arma::vec grad_vals(p), a_vals(p), b_vals(p);
  b_vals.zeros(); a_vals.zeros();

  arma::wall_clock timer;
  timer.tic();

  // Simulate initial times
  arma::uvec inds_off_hp = arma::find(arma::abs(theta) > eps);
  arma::uvec inds_on_hp = arma::find(arma::abs(theta) < eps);
  arma::uvec inds_off_hp_cv = arma::find(arma::abs(cvref) > eps);
  arma::uvec miniu(1);

  // Calcs for ss
  arma::vec C_ss(p);
  for( int i =0; i<p; i++){
    C_ss(i)= n_obs*arma::max(arma::abs(dataX.col(i)));
  }

  // Calcs for cv_ref
  arma::vec en_cv(p), grad_cv(p), C_cv(p);
  int cv_size = inds_off_hp_cv.size();  double ref_rate;
  grad_cv(inds_off_hp_cv) = get_grad_logit(dataX,dataX.cols(inds_off_hp_cv)*cvref(inds_off_hp_cv), datay, inds_off_hp_cv);

  arma::vec X_cv_2(n_obs);
  arma::rowvec dataXi(p);
  for( int i =0; i<n_obs; i++){
    dataXi = dataX.row(i);
    X_cv_2(i) = sqrt(sum(dataXi(inds_off_hp_cv)%dataXi(inds_off_hp_cv)));
  }
  for(unsigned int i =0; i< inds_off_hp_cv.size(); i++){
    C_cv(inds_off_hp_cv(i)) = n_obs*0.25*arma::max(arma::abs(dataX.col(inds_off_hp_cv(i)))%X_cv_2);
  }

  use_cv = check_cv(inds_off_hp_cv, inds_off_hp);
  if(use_cv){
    ref_rate = sum(arma::max(theta(inds_off_hp)%grad_cv(inds_off_hp_cv),arma::zeros(cv_size)));
    a_vals(inds_off_hp) = ref_rate + arma::norm(x(inds_off_hp) - cvref(inds_off_hp),2) * C_cv(inds_off_hp) +
      x(inds_off_hp)%theta(inds_off_hp)/prior_sigma2;
    b_vals(inds_off_hp) = C_cv(inds_off_hp)*sqrt(inds_off_hp.size()) +
      theta(inds_off_hp)%theta(inds_off_hp)/prior_sigma2;// fix later for general
  } else {
    a_vals(inds_off_hp) = C_ss(inds_off_hp) + x(inds_off_hp)%theta(inds_off_hp)/prior_sigma2;
    b_vals(inds_off_hp) = theta(inds_off_hp)%theta(inds_off_hp)/prior_sigma2;
  }

  taus.elem(inds_off_hp) = get_time_off_hp(x.elem(inds_off_hp),theta.elem(inds_off_hp),a_vals.elem(inds_off_hp), b_vals.elem(inds_off_hp));
  taus.elem(inds_on_hp) = zigzag_MP_GaussIID(theta, rj_val, ppi, prior_sigma2);
  while( nEvent < nmax + burn ){
    mini = taus.index_min();
    tau_val = (taus(mini) - t);
    x += tau_val*theta;
    t = taus(mini);
    a_vals += b_vals*tau_val;

    if( std::abs(x(mini)) < eps){
      x(mini) = 0;
      if( std::abs(theta(mini)) < eps ){
        // Exiting the HP
        theta(mini) = 2.0*(R::runif(0,1) < 0.5) - 1;
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
          theta(mini) = 0;
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
      use_cv = check_cv(inds_off_hp_cv, inds_off_hp);
      if(use_cv){
        ref_rate = sum(arma::max(theta(inds_off_hp)%grad_cv(inds_off_hp),arma::zeros(cv_size)));
        a_vals(inds_off_hp) = ref_rate + arma::norm(x(inds_off_hp) - cvref(inds_off_hp),2) * C_cv(inds_off_hp) +
          x(inds_off_hp)%theta(inds_off_hp)/prior_sigma2;
        b_vals(inds_off_hp) = C_cv(inds_off_hp)*sqrt(inds_off_hp.size()) +
          theta(inds_off_hp)%theta(inds_off_hp)/prior_sigma2;// fix later for general
      } else {
        a_vals(inds_off_hp) = C_ss(inds_off_hp) + x(inds_off_hp)%theta(inds_off_hp)/prior_sigma2;
        b_vals(inds_off_hp) = theta(inds_off_hp)%theta(inds_off_hp)/prior_sigma2;
      }

    } else {

      upper = std::max(0.0, a_vals(mini));
      miniu(0) = mini;
      grad_val = get_grad_logit_ss(dataX, x, datay, miniu, cvref, use_cv, grad_cv) +
        x(mini)/prior_sigma2;

      val = theta(mini)*grad_val(0);
      if( R::runif(0,1) < val/upper){
        theta(mini) = -theta(mini);
        if( nEvent >= burn){
          sk_times(nEvent-burn) = t;
          sk_points.col(nEvent-burn) = x;
          sk_theta.col(nEvent-burn) = theta;
        }
        nEvent++;
      }
      if(use_cv){
        ref_rate = sum(arma::max(theta(inds_off_hp)%grad_cv(inds_off_hp_cv),arma::zeros(cv_size)));
        a_vals(inds_off_hp) = ref_rate + arma::norm(x(inds_off_hp) - cvref(inds_off_hp),2) * C_cv(inds_off_hp) +
          x(inds_off_hp)%theta(inds_off_hp)/prior_sigma2;
        b_vals(inds_off_hp) = C_cv(inds_off_hp)*sqrt(inds_off_hp.size()) +
          theta(inds_off_hp)%theta(inds_off_hp)/prior_sigma2;// fix later for general
      } else {
        a_vals(inds_off_hp) = C_ss(inds_off_hp) + x(inds_off_hp)%theta(inds_off_hp)/prior_sigma2;
        b_vals(inds_off_hp) = theta(inds_off_hp)%theta(inds_off_hp)/prior_sigma2;
      }
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

    taus.elem(inds_off_hp) = t + get_time_off_hp(x.elem(inds_off_hp),theta.elem(inds_off_hp),a_vals.elem(inds_off_hp),b_vals.elem(inds_off_hp));
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
