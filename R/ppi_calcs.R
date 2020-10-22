PPI_calc <- function(times, positions, ModelIndices = NULL, MarginalIndicies = NULL, burnin = 1, maxIter=NULL, eps = 1e-12){
  if(!is.null(maxIter)){
    times <- times[1:maxIter]
    positions <- positions[,1:maxIter]
  }
  ## Calc Model probs
  nMod <- length(ModelIndices)
  pipVal <- rep(0,nMod)
  if(!is.null(ModelIndices)){
    for( Mi in 1:nMod ){
      ZeroInds <- which(strsplit(ModelIndices[Mi],'')[[1]] != "1")
      t_dirac = 0
      for(i in burnin:(length(times)-1)){
        if(all(abs(positions[-ZeroInds,i]) > eps) && all(abs(positions[ZeroInds,i]) < eps) &&
           all(positions[ZeroInds,i+1] ==0) && all(abs(positions[-ZeroInds,i+1]) > 0) ){
          t_dirac = t_dirac + times[i+1] - times[i]
        }
      }
      pipVal[Mi] <- t_dirac/(max(times)-times[burnin])
    }
    names(pipVal) <- ModelIndices
  }

  ## Calc Marginal
  nMI <- length(MarginalIndicies)
  margVal <- rep(0, nMI)
  if(!is.null(MarginalIndicies)){
    for( Mi in 1:nMI){
      ZeroInds <- MarginalIndicies[Mi]
      t_dirac <- 0
      for(i in burnin:(length(times)-1)){
        if( positions[ZeroInds,i] == 0 && positions[ZeroInds,i+1] == 0 ){
          t_dirac = t_dirac + times[i+1] - times[i]
        }
      }
      margVal[Mi] <- 1-t_dirac/(max(times)-times[burnin])
    }
    names(margVal) <- MarginalIndicies
  }
  return(list(pip_Models = pipVal, pip_Marginal = margVal))
}

PPI_calc_running <- function(times, positions, MarginalIndicies = NULL, burnin = 1, maxIter=NULL, eps = 1e-12){
  if(!is.null(maxIter)){
    times <- times[1:maxIter]
    positions <- positions[,1:maxIter]
  }
  ## Calc Marginal
  nMI <- length(MarginalIndicies)
  maxIter <- length(times)
  margVal <- matrix(0, nrow = nMI, maxIter -burnin -1)
  if(!is.null(MarginalIndicies)){
    for( Mi in 1:nMI){
      ZeroInds <- MarginalIndicies[Mi]
      t_dirac <- 0
      for(i in (burnin):(maxIter-1)){
        if( positions[ZeroInds,i] == 0 && positions[ZeroInds,i+1] == 0 ){
          t_dirac = t_dirac + times[i+1] - times[i]
        }
        margVal[Mi,i - burnin] <- 1-t_dirac/(times[i+1]-times[burnin])
      }
    }
    #names(margVal) <- MarginalIndicies
  }
  return(list(pip_Marginal = margVal))
}

num_model <- function(times, theta, comptimes, which_seen = F, eps = 1e-10){
  jumpTimes <- which(names(comptimes) != "Switch")
  niter <- length(jumpTimes)
  numvar <- rep(0, niter)
  model_seen <- matrix(abs(theta[,1]) > eps, nrow = 1)
  for( i in 1:niter){
    model_i <- abs(theta[,jumpTimes[i]]) > eps
    if( which_seen ){
      seen <- FALSE
      for(j in 1:nrow(model_seen)){
        if(all(model_seen[j,] == model_i)){
          seen <- TRUE
          break
        }
      }
      if(!seen){
        model_seen <- rbind(model_seen, model_i)
      }
    }
    numvar[i] <- sum(model_i)
  }
  # return number of nonzero coefficients
  return(list(models_seen=model_seen, numer_variables=numvar, comptimes = comptimes[jumpTimes]))
}
num_model_gibbs <- function(gamma, comptimes, which_seen = F, eps = 1e-10){
  niter <- length(gamma[1,])
  numvar <- rep(0, niter)
  model_seen <- matrix(abs(gamma[,1]) > eps, nrow = 1)

  for( i in 1:niter){
    model_i <- abs(gamma[,i]) > eps
    if( which_seen ){
      seen <- FALSE
      for(j in 1:nrow(model_seen)){
        if(all(model_seen[j,] == model_i)){
          seen <- TRUE
          break
        }
      }
      if(!seen){
        model_seen <- rbind(model_seen, model_i)
      }
    }
    numvar[i] <- sum(model_i)
  }
  # return number of nonzero coefficients
  return(list(models_seen=model_seen, numer_variables=numvar,comptimes = comptimes[1:i]))
}


PPI_gibbs_calc <- function(Gamma, ModelIndices = NULL, MarginalIndicies = NULL, burnin = 1, eps = 1e-12){

  ## Calc Model probs
  Gamma <- Gamma[,-c(1:burnin)]
  nMod <- length(ModelIndices)
  pipVal <- rep(0,nMod)
  if(!is.null(ModelIndices)){
    for( Mi in 1:nMod ){
      ZeroInds <- which(strsplit(ModelIndices[Mi],'')[[1]] != "1")
      equalMod <- apply(Gamma, 2, function(ga) all(abs(ga[-ZeroInds])>eps)&all(abs(ga[ZeroInds])<eps))
      pipVal[Mi] <- mean(equalMod)
    }
    names(pipVal) <- ModelIndices
  }
  ## Calc Marginal
  nMI <- length(MarginalIndicies)
  margVal <- rep(0, nMI)
  if(!is.null(MarginalIndicies)){
    for( Mi in 1:nMI){
      margVal[Mi] <- mean(Gamma[MarginalIndicies[Mi],])
    }
    names(margVal) <- MarginalIndicies
  }
  return(list(pip_Models = pipVal, pip_Marginal = margVal))
}
mean_marg <- function(times, positions, thetas, MarginalIndicies = NULL, burnin = 1, eps = 1e-12){
  ## Calc Marginal
  nMI <- length(MarginalIndicies)
  maxIter <- length(times)
  margVal <- matrix(0, nrow = nMI, maxIter -burnin -1)
  pip_Marginal <- rep(0, nMI)
  if(!is.null(MarginalIndicies)){
    for( Mi in 1:nMI){
      mi <- MarginalIndicies[Mi]
      beta_mean <- 0
      total_time <- 0
      for(i in (burnin):(maxIter-1)){
        tauv <- (times[i+1] - times[i])
        total_time <- total_time + tauv
        beta_mean = beta_mean + (tauv*positions[mi,i] + thetas[mi,i]*tauv^2/2)
        margVal[Mi,i - burnin] <- beta_mean/(times[i+1])
        pip_Marginal[Mi] = beta_mean/total_time
      }
    }
  }
  return(list(pip_Marginal_cum = margVal, pip_Marginal=pip_Marginal))
}
cond_marg <- function(times, positions, thetas, theta_c, MarginalIndicies = NULL, burnin = 1, eps = 1e-12){
  ## Calc Marginal
  nMI <- length(MarginalIndicies)
  maxIter <- length(times)
  margVal <- matrix(0, nrow = nMI, maxIter -burnin -1)
  pip_Marginal <- rep(0, nMI)
  if(!is.null(MarginalIndicies)){
    for( Mi in 1:nMI){
      mi <- MarginalIndicies[Mi]
      beta_mean <- 0
      total_time <- 0
      for(i in (burnin):(maxIter-1)){
        tauv <- (times[i+1] - times[i])
        if(all(abs(theta_c-abs(thetas[,i])) <1e-10)){
          total_time <- total_time + tauv
          beta_mean = beta_mean + (tauv*positions[mi,i] + thetas[mi,i]*tauv^2/2)
          margVal[Mi,i - burnin] <- beta_mean/(times[i+1])
          pip_Marginal[Mi] = beta_mean/total_time
        }
      }
    }
  }
  return(list(pip_Marginal_cum = margVal, pip_Marginal=pip_Marginal))
}
post_cummean <- function(times, positions, thetas, MarginalIndicies = NULL, burnin = 1, maxIter=NULL, eps = 1e-12){
  if(!is.null(maxIter)){
    times <- times[burnin:maxIter] - times[burnin]
    positions <- positions[,burnin:maxIter]
    thetas <- thetas[,burnin:maxIter]
  }
  ## Calc Marginal
  nMI <- length(MarginalIndicies)
  maxIter <- length(times)
  margVal <- matrix(0, nrow = nMI, maxIter -burnin -1)
  if(!is.null(MarginalIndicies)){
    for( Mi in 1:nMI){
      mi <- MarginalIndicies[Mi]
      beta_mean <- 0
      for(i in (burnin):(maxIter-1)){
        tauv <- (times[i+1] - times[i])
        beta_mean = beta_mean + tauv*positions[mi,i] + thetas[mi,i]*tauv^2/2
        margVal[Mi,i - burnin] <- beta_mean/(times[i+1])
      }
    }
    #names(margVal) <- MarginalIndicies
  }
  return(list(pip_Marginal = margVal))
}
