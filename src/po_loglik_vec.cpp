#include <Rcpp.h>
using namespace Rcpp;

//' Per-row proportional-odds log-likelihoods
//'
//' Computes per-observation log-likelihood contributions for a
//' proportional-odds model (ac = FALSE).
//'
//' @param theta Numeric vector of parameters:
//'   \eqn{(\alpha_1, \ldots, \alpha_{K-1}, \beta_1, \ldots, \beta_p)}.
//' @param Xmat Numeric matrix of predictors (n x p).
//' @param cumY Numeric matrix of cumulative indicators (n x K).
//' @param eps Small constant for numerical stability.
//' @return Numeric vector of length n with per-row log-likelihoods.
// [[Rcpp::export]]
NumericVector po_loglik_vec_cpp(const NumericVector& theta,
                                const NumericMatrix& Xmat,
                                const NumericMatrix& cumY,
                                const double eps = 1e-12) {
  int n = Xmat.nrow();
  int p = Xmat.ncol();
  int k = cumY.ncol();
  if (k < 2) stop("cumY must have at least 2 columns.");
  
  int k1 = theta.size() - p; // #cutpoints
  if (k1 != (k - 1)) stop("length(theta) - ncol(Xmat) must equal k-1.");
  
  // Split theta
  std::vector<double> alpha(k1);
  for (int r = 0; r < k1; ++r) alpha[r] = theta[r];
  std::vector<double> beta(p);
  for (int j = 0; j < p; ++j) beta[j] = theta[k1 + j];
  
  NumericVector out(n);
  
  for (int i = 0; i < n; ++i) {
    // linear predictor: X_i %*% beta
    double eta_lin = 0.0;
    for (int j = 0; j < p; ++j) eta_lin += Xmat(i, j) * beta[j];
    
    double loglik = 0.0;
    for (int r = 0; r < k1; ++r) {
      double eta = alpha[r] + eta_lin;
      // logistic
      double gamma = 1.0 / (1.0 + std::exp(-eta));
      if (gamma < eps) gamma = eps;
      if (gamma > 1.0 - eps) gamma = 1.0 - eps;
      
      // gamma_{r+1}
      double gamma_next;
      if (r < k1 - 1) {
        double eta_next = alpha[r + 1] + eta_lin;
        gamma_next = 1.0 / (1.0 + std::exp(-eta_next));
        if (gamma_next < eps) gamma_next = eps;
        if (gamma_next > 1.0 - eps) gamma_next = 1.0 - eps;
      } else {
        gamma_next = 1.0; // last one
      }
      
      double ph  = std::log(gamma) - std::log(std::max(gamma_next - gamma, eps));
      double gph = std::log(1.0 + std::exp(ph));
      
      double Yj   = cumY(i, r); // cumYi[r]
      double Yjp1 = (r < k1 - 1) ? cumY(i, r + 1) : 1.0;
      
      loglik += Yj * ph - Yjp1 * gph;
    }
    out[i] = loglik;
  }
  
  return out;
}
