#' Generate samples for PDMP trajectory
#'
#' Get samples from PDMP trajectories taking a fixed time discretisation.
#'
#' @param positions Matrix of positions from the PDMP trajectory, each column should correspond to a position
#' @param times Vector of event times from the PDMP trajectory
#' @param nsample Number of desired samples from the PDMP trajectory
#' @param theta Optional Matrix of velocities from the PDMP trjectory, each column should correspond to a velocity
#' @return Returns a list with the following objects:
#' @return \code{x}: Matrix of extracted samples of the position (x) taken using a fixed time discretisation of the PDMP
#' @return \code{theta}: Matrix of extracted samples of the velocity (theta) taken using a fixed time discretisation of the PDMP
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
#' samples <- gen_sample(zigzag_fit$positions, zigzag_fit$times, 10^4)
#'
#' plot(zigzag_fit$positions[1,],zigzag_fit$positions[2,], type = 'l', xlab = 'x1', ylab = 'x2')
#' points(samples$xx[1,], samples$xx[2,], col='red', pch=20)
#'
gen_sample <- function(positions, times, nsample, theta=NULL, burn = 1){
  if(is.null(dim(positions))) positions <- matrix(positions, nrow = 1)
  positions <- positions[,burn:length(times)]
  times <- times[burn:length(times)] - times[burn]
  nsteps <- length(times)
  Tmax <- times[nsteps]
  dt <- Tmax/(nsample + 1)
  t = dt
  t0 = times[1]
  x0 = positions[,1]
  thetas<-matrix(0, nrow = length(x0), ncol = nsample)
  if(!is.null(theta)){
    theta0 <-theta[,1]
  }
  samples <- matrix(0, nrow = length(x0), ncol = nsample)
  n <- 0

  for(i in 2:nsteps){
    x1 = positions[,i]
    t1 = times[i]
    if(!is.null(theta)) {theta1 = theta[,i]}
    while(t < t1 && n < nsample){
      n <- n+1
      samples[,n] <- x0 + (x1 - x0)*(t-t0)/(t1-t0)
      if(!is.null(theta)) {thetas[,n] <- theta0}
      t <- t + dt
    }
    x0 = x1; t0 = t1; theta0 = if(!is.null(theta)) {theta1}
  }
  return(list(x = samples, theta = thetas))
}
