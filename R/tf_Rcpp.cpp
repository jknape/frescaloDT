#include <Rcpp.h>
using namespace Rcpp;
#include <Rcpp.h>

/*
while(k <= kmax & !conv) {
  estval = 1 - exp(-plog * tf)
  esttot = sum(wgt * estval)
  estvar = sum(wgt^2 * estval * (1-estval))
  if (abs(sptot-esttot) < 0.0005) {
    conv = TRUE
  } else {
    tf=tf*sptot/(esttot+0.0000001)
  }
  k = k +1
}

 Need: tf, esttot, estvar, k
*/


// [[Rcpp::export]]
NumericVector tf_iter(
             double tf,
             const NumericVector wgt,
             const NumericVector plog,
             double sptot,
             int kmax = 1000,
             double tol = 0.0005) {
  // NOTE: Rcpp doesn't allow default reference arguments like in plain C++.
  // Caller must provide 'conv' bool variable.
 // NumericVector out (3);
  const int n = wgt.size();
  const double eps = 1e-7;
  bool conv = false;
  double one_minus_exp;
  double w;
  int k = 1;
  double esttot = 0.0;
  double estvar = 0.0;
  while (k <= kmax && !conv) {
    esttot = 0.0;
    estvar = 0.0;
    for (int i = 0; i < n; ++i) {
      one_minus_exp = 1.0 - exp(-plog[i] * tf);
      w = wgt[i];
      esttot += w * one_minus_exp;
      estvar += (w * w) * one_minus_exp * (1.0 - one_minus_exp);
    }
    if (std::abs(sptot - esttot) < tol) {
      conv = true;
  } else {
      tf = tf * sptot / (esttot + eps);
    ++k;
  }
  }
 //return List::create(tf,esttot, k, estvar);
 return {tf, esttot, estvar,(double) k};
  //return out;
  // If caller needs esttot/estvar/estval, they should provide containers and modify this function accordingly.
}
