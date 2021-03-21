#' Plot PDMP dynamics and samples for posterior distributions
#'
#' Plot marginal densities and joint pairs plots for trajectories and samples of PDMP samplers and optionally MCMC samples for comparison.
#' Care should be taken when interpreting marginal KDE estimates on the diagonal as the bandwidth of the KDE has an impact on how the Dirac spike is visualised.
#' @param pdmp_res List of positions, times and velocities returned from a PDMP sampler
#' @param coords Vector of coordinates to plot the marginal and joint distributions
#' @param inds Vector of indices of the PDMP trajectories to plot.
#' @param nsamples Number of samples to generate and use for marginal density estimates of the PDMP methods
#' @param burn Percentage of events to use as burn-in. Should be between 0 and 1.
#' @param mcmc_samples Optional Matrix of samples from an MCMC method. Each row should be a sample.
#' @param pch The graphics parameter for off diagonal plots. Default is 20.
#' @param cols Colours to be used for plotting the PDMPs and MCMC samples (in order).
#' @return Generates a plot of the marginal density on the diagonal and pairs plots of the trajectories
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
#' gibbs_fit <- gibbs_logit(maxTime = 1, dataX = data$dataX, datay = data$dataY,
#'                          prior_sigma2 = 10,beta = rep(0,p), gamma = rep(0,p),
#'                          ppi = ppi)
#'\dontrun{
#' plot_pdmp(zigzag_fit, coords = 1:2, inds = 1:10^3,burn = .1,
#'           nsamples = 1e4, mcmc_samples = t(gibbs_fit$beta*gibbs_fit$gamma))
#'}
#'
plot_pdmp <- function(pdmp_res, coords = 1:2, inds = 1:10^3, nsamples = 10^3,
                      burn = 0.1, mcmc_samples=NULL, pch = 20, cols = NULL){
  ndim <- length(coords)
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mfrow = c(ndim,ndim),
      mar = rep(2,4))
  samples <- gen_sample(pdmp_res$positions, pdmp_res$times,theta = pdmp_res$theta,
                        nsample = nsamples, burn = burn*length(pdmp_res$times))

  plot_mcmc <- !is.null(mcmc_samples)
  leg_names <- if(plot_mcmc) c("PDMP", "MCMC") else c("PDMP")
  colres <-  if(is.null(cols)){1:(1+plot_mcmc)} else {cols}

  for( i in 1:ndim ){
    for(j in 1:ndim ){
      if(i == j & nsamples > 0){
        ds <- density(samples$x[coords[i],])
        plot(ds, main='',xlab='', ylab='', col = colres[1])
        if(plot_mcmc){
          lines(density(mcmc_samples[,coords[i]], bw = ds$bw), col = colres[2])
        }
      }
      if( i != j ){
        xrange <- range(c(samples$x[coords[i],], pdmp_res$positions[coords[i],inds]))
        yrange <- range(c(samples$x[coords[j],], pdmp_res$positions[coords[j],inds]))
        plot(pdmp_res$positions[coords[i],inds], pdmp_res$positions[coords[j],inds],
             xlim = xrange, ylim = yrange, type = 'l', main = '')
        points(samples$x[coords[i],], samples$x[coords[j],], col = colres[1], pch = pch)
        if(plot_mcmc){
          points(mcmc_samples[,coords[i]], mcmc_samples[,coords[j]], col = colres[2], pch = pch)
        }
      }
    }
  }

  legend('topleft',legend = leg_names, col = colres, lwd = 1 )
}

#' Plot multiple PDMP dynamics and MCMC samples for posterior distributions
#'
#' Plots to compare PDMP samplers and optionally MCMC samples.
#' @param list_pdmp List of PDMP sampler trajectories to plot
#' @param coords Vector of coordinates to plot the marginal and joint distributions
#' @param inds Vector of indices of the PDMP trajectories to plot.
#' @param nsamples Number of samples to generate and use for marginal density estimates of the PDMP methods
#' @param burn Percentage of events to use as burn-in. Should be between 0 and 1, default 0.1.
#' @param mcmc_samples Optional Matrix of samples from an MCMC method. Each row should be a sample.
#' @param pch The graphics parameter for off diagonal plots. Default is 20.
#' @param cols Colours to be used for plotting the PDMPs and MCMC samples (in order).
#' @return Generates a plot of the marginal density on the diagonal and pairs plots of the trajectories
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
#' zigzag_fit <- zigzag_logit(maxTime = 1, dataX = data$dataX,datay =data$dataY,
#'                            prior_sigma2 = 10,theta0 = rep(0, p),
#'                            x0 = rep(0, p), rj_val = 0.6,
#'                            ppi = ppi)
#'\dontrun{
#' bps_fit <- bps_n_logit(maxTime = 1, dataX = data$dataX, datay = data$dataY,
#'                        prior_sigma2 = 10, theta0 = rep(0, p), x0 = rep(0, p),
#'                        ref = 0.1, rj_val = 0.6,ppi = ppi)
#'
#' plot_pdmp_multiple(list(zz=zigzag_fit,bps=bps_fit), coords = 1:2, inds = 1:10^3,
#'                       nsamples = 1e4, burn = .1)
#'}
#'
plot_pdmp_multiple <- function(list_pdmp, coords = 1:2, inds = 1:10^3,
                               nsamples = 10^3, burn = 0.1, mcmc_samples=NULL, pch = 20, cols = NULL){
  ndim <- length(coords)
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mfrow = c(ndim,ndim),
      mar = rep(2,4))
  nres <- length(list_pdmp)
  if(is.null(names(list_pdmp))){
    names(list_pdmp) <- 1:nres
  }
  samples <- lapply(list_pdmp, function(res)
    rjpdmp::gen_sample(res$positions, res$times,theta = res$theta,
                       nsample = nsamples, burn = burn*length(res$times)))

  plot_mcmc <- !is.null(mcmc_samples)
  leg_names <- if(plot_mcmc) c(names(list_pdmp), "MCMC") else names(list_pdmp)

  colres <-  if(is.null(cols)){1:(nres+plot_mcmc)} else {cols}
  for( i in 1:ndim ){
    for(j in 1:ndim ){
      if(i == j & nsamples > 0){
        ds <- density(samples[[1]]$x[coords[i],])
        plot(ds, main='',xlab='', ylab='', col = colres[1])
        for( r in 2:nres){
          lines(density(samples[[r]]$x[coords[i],], bw = ds$bw), col = colres[r])
        }
        if(plot_mcmc){
          lines(density(mcmc_samples[,coords[i]]), col = colres[nres+plot_mcmc])
        }
      }
      if( i != j ){
        xrange <- range(sapply(1:nres, function(rs) range(c(samples[[rs]]$x[coords[i],], list_pdmp[[rs]]$positions[coords[i],inds]))))
        yrange <- range(sapply(1:nres, function(rs) range(c(samples[[rs]]$x[coords[j],], list_pdmp[[rs]]$positions[coords[j],inds]))))

        plot(list_pdmp[[1]]$positions[coords[i],inds], list_pdmp[[1]]$positions[coords[j],inds],xlim = xrange,
             ylim = yrange, type = 'l', main = '', col = colres[1])
        points(samples[[1]]$x[coords[i],], samples[[1]]$x[coords[j],], col = colres[1], pch = pch)
        for( r in 2:nres){
          lines(list_pdmp[[r]]$positions[coords[i],inds], list_pdmp[[r]]$positions[coords[j],inds], col = colres[r])
          points(samples[[r]]$x[coords[i],], samples[[r]]$x[coords[j],], col = colres[r], pch = pch)
        }

        if(plot_mcmc){
          points(mcmc_samples[,coords[i]], mcmc_samples[,coords[j]], col = colres[nres+plot_mcmc], pch = pch)
        }
      }
    }
  }
  legend('topleft',legend = leg_names, col = colres, lwd = 1 )
}
