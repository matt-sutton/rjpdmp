#ifndef GET_INDEX_H
#define GET_INDEX_H

// [[Rcpp::export]]
arma::uvec get_index(int nsamp, int n_obs){
  arma::uvec indices(nsamp);
  for( int i = 0; i < nsamp; i++){
    indices(i) = std::floor(R::runif(0,1)*n_obs);
  }
  return(indices);
}

#endif
