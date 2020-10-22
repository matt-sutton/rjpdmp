

plot_pdmp <- function(pdmp_res, margins = 1:2, inds = 1:10^3, nsamples = 10^3, burn = 0.1, mcmc_samples=NULL){
  ndim <- length(margins)
  par(mfrow = c(ndim,ndim),
      mar = rep(2,4))
  samples <- gen_sample(pdmp_res$positions, pdmp_res$times,theta = pdmp_res$theta,
                                   nsample = nsamples, burn = burn*length(pdmp_res$times))
  for( i in 1:ndim ){
    for(j in 1:ndim ){
      if(i == j & nsamples > 0){
        plot(density(samples$xx[margins[i],]), main='',xlab='', ylab='')
        if(!is.null(mcmc_samples)){
          lines(density(mcmc_samples[,margins[i]]), col = 'blue')
        }
      }
      if( i != j ){
        xrange <- range(c(samples$xx[margins[i],], pdmp_res$positions[margins[i],inds]))
        yrange <- range(c(samples$xx[margins[j],], pdmp_res$positions[margins[j],inds]))
        plot(pdmp_res$positions[margins[i],inds], pdmp_res$positions[margins[j],inds],xlim = xrange, ylim = yrange, type = 'l', main = '')
        points(samples$xx[margins[i],], samples$xx[margins[j],], col = 'red', pch = '.')
        if(!is.null(mcmc_samples)){
          points(mcmc_samples[,margins[i]], mcmc_samples[,margins[j]], col = 'blue', pch = '.')
        }
      }
    }
  }
}

plot_pdmp_multiple <- function(list_pdmp, margins = 1:2, inds = 1:10^3,
         nsamples = 10^3, burn = 0.1, mcmc_samples=NULL){
  ndim <- length(margins)
  par(mfrow = c(ndim,ndim),
      mar = rep(2,4))
  nres <- length(list_pdmp)
  samples <- lapply(list_pdmp, function(res)
    gen_sample(res$positions, res$times,theta = res$theta,
                          nsample = nsamples, burn = burn*length(res$times)))
  colres <- rainbow(nres)
  for( i in 1:ndim ){
    for(j in 1:ndim ){
      if(i == j & nsamples > 0){
        plot(density(samples[[1]]$xx[margins[i],]), main='',xlab='', ylab='', col = colres[1])
        for( r in 2:nres){
          lines(density(samples[[r]]$xx[margins[i],]), col = colres[r])
        }
        if(!is.null(mcmc_samples)){
          lines(density(mcmc_samples[,margins[i]]), col = 'blue')
        }
      }
      if( i != j ){
        xrange <- range(sapply(1:nres, function(rs) range(c(samples[[rs]]$xx[margins[i],], list_pdmp[[rs]]$positions[margins[i],inds]))))
        yrange <- range(sapply(1:nres, function(rs) range(c(samples[[rs]]$xx[margins[j],], list_pdmp[[rs]]$positions[margins[j],inds]))))

        plot(list_pdmp[[1]]$positions[margins[i],inds], list_pdmp[[1]]$positions[margins[j],inds],xlim = xrange,
             ylim = yrange, type = 'l', main = '', col = colres[1])
        points(samples[[1]]$xx[margins[i],], samples[[1]]$xx[margins[j],], col = colres[1], pch = '.')
        for( r in 2:nres){
          lines(list_pdmp[[r]]$positions[margins[i],inds], list_pdmp[[r]]$positions[margins[j],inds], col = colres[r])
          points(samples[[r]]$xx[margins[i],], samples[[r]]$xx[margins[j],], col = colres[r], pch = '.')
        }

        if(!is.null(mcmc_samples)){
          points(mcmc_samples[,margins[i]], mcmc_samples[,margins[j]], col = 'blue', pch = '.')
        }
      }
    }
  }
}
