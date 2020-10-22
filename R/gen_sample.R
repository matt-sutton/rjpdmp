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
  return(list(xx = samples, theta = thetas))
}
