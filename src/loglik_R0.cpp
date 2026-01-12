#include <Rcpp.h>
using namespace Rcpp;

// Fast R0 loglik aggregation: sum_{i} log( sum_{k} prob_po[i,k] * sum_j B_ij p_kj )
// [[Rcpp::export]]
double loglik_R0_cpp(const NumericVector& prob_po,
                     const NumericMatrix& B_mat,
                     const NumericMatrix& p_mat,
                     const IntegerVector& id_vec, // 1..n_R0
                     const int n_R0) {
  const int N = B_mat.nrow();
  const int S = B_mat.ncol();
  NumericVector lik_id(n_R0);
  
  for (int i = 0; i < N; ++i) {
    double row_sum = 0.0;
    const double w = prob_po[i];
    for (int s = 0; s < S; ++s) row_sum += B_mat(i, s) * p_mat(i, s);
    const int id0 = id_vec[i] - 1;
    lik_id[id0] += w * row_sum;
  }
  
  double out = 0.0;
  for (int id0 = 0; id0 < n_R0; ++id0) {
    const double v = (lik_id[id0] > 0.0) ? lik_id[id0] : 1.0;
    out += std::log(v);
  }
  return out;
}
