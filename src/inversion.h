#ifndef INVERSION_H
#define INVERSION_H

// Simulate from a + bs
// [[Rcpp::export]]
double linear_inv_t(double a, double b, double u) {
  if( a < 0 & b > 0){
    return(-a/b+sqrt(-2*log(u)/b));

  } else if ( a > 0 & b < 0 ){
    if(-pow(a,2)/b + pow(-a,2)/(2*b) >= -log(u)){
      return(-a/b -sqrt(pow(a/b,2) - 2*log(u)/b));
    } else {
      return(R_PosInf);
    }

  } else if ( a >= 0 & b > 0 ){
    return(-a/b + sqrt(pow(a/b,2) - 2*log(u)/b));

  } else if ( a > 0 & b == 0 ){
    return(-log(u)/a);

  } else if ( a <= 0 & b <= 0 ){
    return(R_PosInf);
  }
}

// [[Rcpp::export]]
double get_time_off_hp_ind(double x, double theta,
                          double a, double b) {
  double hit_time = R_PosInf;
  if( (x*theta < 1e-14) && (std::abs(x) > 1e-14)){
    hit_time = -x/theta;
  }
  return(std::min(hit_time,linear_inv_t(a,b,R::runif(0,1))));
}

// [[Rcpp::export]]
arma::vec get_time_off_hp(arma::vec x, arma::vec theta,
                          arma::vec a, arma::vec b) {
  int p = x.size();

  arma::vec hit_times = R_PosInf*arma::ones(p,1);
  arma::uvec nzero = find((x%theta < 1e-14) && (abs(x) > 1e-14));
  hit_times.elem(nzero) = -x.elem(nzero)/theta.elem(nzero);

  arma::vec times(p);
  for( int i = 0; i < p; i++ ){
    times[i] = std::min(hit_times[i],linear_inv_t(a[i],b[i],R::runif(0,1)));
  }
  return(times);
}

// [[Rcpp::export]]
arma::vec get_hit_times(arma::vec x, arma::vec theta) {
  int p = x.size();
  arma::vec hit_times = R_PosInf*arma::ones(p,1);
  arma::uvec nzero = find((x%theta < 1e-14) && (abs(x) > 1e-14));
  hit_times.elem(nzero) = -x.elem(nzero)/theta.elem(nzero);
  return(hit_times);
}

#endif
