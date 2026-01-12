#include <Rcpp.h>
using namespace Rcpp;

// E-step fused kernel
//    Inputs (length N = n_R0 * m):
//     - B_mat: N x S, p_mat: N x S, prob_po: N
//     - id_vec: 1..n_R0 (group for denominator), k_vec: 1..m (group for p update)
//    Returns:
//     - q_vec: length N
//     - p_kj_num_R0: m x S (rowsum of normalized psi by k)
// [[Rcpp::export]]
List estep_qvec_cpp(const NumericMatrix& B_mat,
                    const NumericMatrix& p_mat,
                    const NumericVector& prob_po,
                    const IntegerVector& id_vec, // 1..n_R0
                    const IntegerVector& k_vec,  // 1..m
                    const int n_R0,
                    const int m) {
  const int N = B_mat.nrow();
  const int S = B_mat.ncol();
  
  NumericVector den_id(n_R0);
  // first pass: accumulate denominator per id
  for (int i = 0; i < N; ++i) {
    const double w = prob_po[i];
    const int id0 = id_vec[i] - 1; // zero-based
    double sum_row = 0.0;
    for (int s = 0; s < S; ++s) sum_row += B_mat(i, s) * p_mat(i, s) * w;
    den_id[id0] += sum_row;
  }
  for (int id0 = 0; id0 < n_R0; ++id0) if (den_id[id0] <= 0.0) den_id[id0] = 1.0;
  
  NumericVector q_vec(N);
  NumericMatrix pkj_num_R0(m, S);
  
  // second pass: normalized psi -> q and rowsum by k
  for (int i = 0; i < N; ++i) {
    const double w = prob_po[i];
    const int id0 = id_vec[i] - 1;
    const int k0  = k_vec[i] - 1;
    const double denom = den_id[id0];
    
    double q_sum = 0.0;
    for (int s = 0; s < S; ++s) {
      const double psi_norm = (B_mat(i, s) * p_mat(i, s) * w) / denom;
      q_sum += psi_norm;
      pkj_num_R0(k0, s) += psi_norm;
    }
    q_vec[i] = q_sum;
  }
  
  return List::create(_["q_vec"] = q_vec, _["p_kj_num_R0"] = pkj_num_R0);
}
