#' Calculate posterior pobabilities of inclusion based on PDMP trajectories
#'
#' Calculate either marginal probabilities of inclusions or posterior probabilities of specific models.
#'
#' @param times Vector of event times from the PDMP trajectory
#' @param positions Matrix of positions from the PDMP trajectory, each column should correspond to a position
#' @param models Optional Matrix of indicies where rows correspond to models. Will return proabilities of each model \code{prob_mod}.
#' @param marginals Optional Vector of indices to calculate the marginal probabilities of inclusion. Will return proabilities of inclusion for variable index \code{marginal_prob}.
#' @param burnin Number of events to use as burnin
#' @return Returns a list with the following objects:
#' @return \code{prob_mod}: Vector of posterior model probabilities based on the PDMP trajectories
#' @return \code{marginal_prob}: Vector of marginal probabilities for inclusion
#' @examples
#' generate.logistic.data <- function(beta, n.obs, Sig) {
#' p <- length(beta)
#' dataX <- MASS::mvrnorm(n=n.obs,mu=rep(0,p),Sigma=Sig)
#' vals <- dataX %*% as.vector(beta)
#' generateY <- function(p) { rbinom(1, 1, p)}
#' dataY <- sapply(1/(1 + exp(-vals)), generateY)
#' return(list(dataX = dataX, dataY = dataY))
#' }
#'
#' n <- 15
#' p <- 25
#' beta <- c(1, rep(0, p-1))
#' Siginv <- diag(1,p,p)
#' Siginv[1,2] <- Siginv[2,1] <- 0.9
#' set.seed(1)
#' data <- generate.logistic.data(beta, n, solve(Siginv))
#' ppi <- 2/p
#'
#' zigzag_fit <- zigzag_logit(maxTime = 1, dataX = data$dataX, datay = data$dataY,
#'                            prior_sigma2 = 10,theta0 = rep(0, p), x0 = rep(0, p), rj_val = 0.6,
#'                            ppi = ppi)
#' \dontrun{
#' a <- models_visited(zigzag_fit$theta)
#'
#' # Work out probability of top 10 most visited models and all marginal inclusion probabilities
#' # specific model probabilities become trivially small for large dimensions
#' b <- model_probabilities(zigzag_fit$times, zigzag_fit$positions,
#'                          models = a[1:10,1:p], marginals=1:p)
#' }
#'
#'
model_probabilities <- function(times, positions, models = NULL, marginals = NULL, burnin = 1){
  eps = 1e-12
  ## Calculate model probabilities
  prob_mod <- NULL
  if(!is.null(models)){
    if(is.null(nrow(models))){
      models <- matrix(models, nrow = 1)
    }
    nMod <- nrow(models)
    prob_mod <- rep(0,nMod)
    for( Mi in 1:nMod ){
      t_in_model = 0
      in_model <- FALSE
      t_at_model <- 0
      for(i in burnin:(length(times)-1)){
        if(all(abs(abs(thetas[,i]) - abs(models[Mi,])) < eps) & !(in_model)){
          t_at_model <- times[i]
          in_model <- TRUE
        }
        if(!all(abs(abs(thetas[,i]) - abs(models[Mi,])) < eps) & (in_model)){
          t_in_model <- t_in_model + (times[i] - t_at_model)
          in_model <- FALSE
        }
      }
      if(in_model){
        t_in_model <- t_in_model + (times[i] - t_at_model)
      }
      prob_mod[Mi] <- exp(log(t_in_model) - log(max(times)-times[burnin]))
    }
    names(prob_mod) <- apply(models,1,function(s) paste(s, collapse = ''))
  }

  ## Calculate Marginal probabilities
  nMI <- length(marginals)
  marginal_prob_inclusion <- rep(0, nMI)
  if(!is.null(marginals)){
    for( Mi in 1:nMI){
      ZeroInds <- marginals[Mi]
      t_dirac <- 0
      for(i in burnin:(length(times)-1)){
        if( thetas[ZeroInds,i] == 0 && thetas[ZeroInds,i+1] == 0 ){
          t_dirac = t_dirac + times[i+1] - times[i]
        }
      }
      marginal_prob_inclusion[Mi] <- 1-t_dirac/(max(times)-times[burnin])
    }
    names(marginal_prob_inclusion) <- marginals
  }
  return(list(prob_mod = sort(prob_mod,decreasing = T), marginal_prob = marginal_prob_inclusion))
}

#' Count the number of times a model is visited
#'
#' @param thetas Vector of model indicies from the PDMP trajectory or samples from an MCMC sampler
#' @return Returns a Matrix with rows corresponding to models and a final column corresponding to the number of times the model is visited
#' @examples
#' generate.logistic.data <- function(beta, n.obs, Sig) {
#' p <- length(beta)
#' dataX <- MASS::mvrnorm(n=n.obs,mu=rep(0,p),Sigma=Sig)
#' vals <- dataX %*% as.vector(beta)
#' generateY <- function(p) { rbinom(1, 1, p)}
#' dataY <- sapply(1/(1 + exp(-vals)), generateY)
#' return(list(dataX = dataX, dataY = dataY))
#' }
#'
#' n <- 15
#' p <- 25
#' beta <- c(1, rep(0, p-1))
#' Siginv <- diag(1,p,p)
#' Siginv[1,2] <- Siginv[2,1] <- 0.9
#' set.seed(1)
#' data <- generate.logistic.data(beta, n, solve(Siginv))
#' ppi <- 2/p
#'
#' zigzag_fit <- zigzag_logit(maxTime = 1, dataX = data$dataX, datay = data$dataY,
#'                            prior_sigma2 = 10,theta0 = rep(0, p), x0 = rep(0, p),
#'                            rj_val = 0.6, ppi = ppi)
#' \dontrun{
#' models_visited(zigzag_fit$theta)
#' }
#'
models_visited <- function(thetas){
  . <- NULL # necessity for data.table package dependency
  df <- data.table::data.table(t(abs(thetas)) > 1e-10)
  df <- as.matrix(df[,.(COUNT = .N), by = names(df)])
  df <- df[order(df[,NCOL(df)], decreasing = T),]
  vnames <- 1:(NCOL(df)-1)
  colnames(df) <- c(vnames, 'count')
  return(df)
}

#' Calculate the marginal mean
#'
#' @param times Vector of event times from the PDMP trajectory
#' @param positions Matrix of positions from the PDMP trajectory, each column should correspond to a position
#' @param thetas Matrix of PDMP velocities
#' @param marginals Vector of indices to calculate the marginal means.
#' @param burnin Number of events to use as burnin
#' @return Returns the posterior mean of the parameter estimated using the PDMP trajectories.
#' @examples
#' generate.logistic.data <- function(beta, n.obs, Sig) {
#' p <- length(beta)
#' dataX <- MASS::mvrnorm(n=n.obs,mu=rep(0,p),Sigma=Sig)
#' vals <- dataX %*% as.vector(beta)
#' generateY <- function(p) { rbinom(1, 1, p)}
#' dataY <- sapply(1/(1 + exp(-vals)), generateY)
#' return(list(dataX = dataX, dataY = dataY))
#' }
#'
#' n <- 15
#' p <- 25
#' beta <- c(1, rep(0, p-1))
#' Siginv <- diag(1,p,p)
#' Siginv[1,2] <- Siginv[2,1] <- 0.9
#' set.seed(1)
#' data <- generate.logistic.data(beta, n, solve(Siginv))
#' ppi <- 2/p
#'
#' zigzag_fit <- zigzag_logit(maxTime = 1, dataX = data$dataX, datay = data$dataY,
#'                            prior_sigma2 = 10,theta0 = rep(0, p), x0 = rep(0, p), rj_val = 0.6,
#'                            ppi = ppi)
#' \dontrun{
#' b <- marginal_mean(zigzag_fit$times, zigzag_fit$positions, zigzag_fit$theta, marginals=1:p)
#' }
#'
marginal_mean <- function(times, positions, thetas, marginals = NULL, burnin = 1){
  ## Calc Marginal
  eps = 1e-12
  nMI <- length(marginals)
  maxIter <- length(times)
  marg_mean <- rep(0, nMI)
  if(!is.null(marginals)){
    for( Mi in 1:nMI){
      mi <- marginals[Mi]
      beta_mean <- 0
      total_time <- 0
      for(i in (burnin):(maxIter-1)){
        tauv <- (times[i+1] - times[i])
        total_time <- total_time + tauv
        beta_mean = beta_mean + (tauv*positions[mi,i] + thetas[mi,i]*tauv^2/2)
        marg_mean[Mi] = beta_mean/total_time
      }
    }
  }
  return(marg_mean)
}

#' Calculate the mean conditioned on being in a specific model
#'
#' @param times Vector of event times from the PDMP trajectory
#' @param positions Matrix of positions from the PDMP trajectory, each column should correspond to a position
#' @param thetas Matrix of PDMP velocities
#' @param theta_c Vector indicating the model to condition on, 1s for active variables and zeros for inactive variables
#' @param burnin Number of events to use as burnin
#' @return Returns the mean conditioned on being in model theta_c estimated using the PDMP trajectories.
#' @examples
#' generate.logistic.data <- function(beta, n.obs, Sig) {
#' p <- length(beta)
#' dataX <- MASS::mvrnorm(n=n.obs,mu=rep(0,p),Sigma=Sig)
#' vals <- dataX %*% as.vector(beta)
#' generateY <- function(p) { rbinom(1, 1, p)}
#' dataY <- sapply(1/(1 + exp(-vals)), generateY)
#' return(list(dataX = dataX, dataY = dataY))
#' }
#'
#' n <- 15
#' p <- 25
#' beta <- c(1, rep(0, p-1))
#' Siginv <- diag(1,p,p)
#' Siginv[1,2] <- Siginv[2,1] <- 0.9
#' set.seed(1)
#' data <- generate.logistic.data(beta, n, solve(Siginv))
#' ppi <- 2/p
#'
#' zigzag_fit <- zigzag_logit(maxTime = 1, dataX = data$dataX, datay = data$dataY,
#'                            prior_sigma2 = 10,theta0 = rep(0, p), x0 = rep(0, p), rj_val = 0.6,
#'                            ppi = ppi)
#' \dontrun{
#' b <- cond_mean(zigzag_fit$times, zigzag_fit$positions, zigzag_fit$theta, theta_c = c(1,rep(0,p-1)))
#' }
#'
cond_mean <- function(times, positions, thetas, theta_c, burnin = 1){
  eps = 1e-12
  nonZeroindicies <- which(abs(theta_c) >1e-10)
  nZ <- length(nonZeroindicies)
  maxIter <- length(times)
  cond_mean <- rep(0, nZ)
  for( Zi in 1:nZ){
    zi <- nonZeroindicies[Zi]
    beta_mean <- 0
    total_time <- 0
    for(i in (burnin):(maxIter-1)){
      tauv <- (times[i+1] - times[i])
      if(all(abs(theta_c-abs(thetas[,i])) <1e-10)){
        total_time <- total_time + tauv
        beta_mean = beta_mean + (tauv*positions[zi,i] + thetas[zi,i]*tauv^2/2)
        cond_mean[Zi] = beta_mean/total_time
      }
    }
  }
  return(cond_mean)
}
