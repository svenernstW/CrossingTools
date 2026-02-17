// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include "parallel_backend.h"  // portable threading: OpenMP if present, else RcppParallel

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP cpp_u_from_from_g_simple(const NumericMatrix& M,
                       const NumericMatrix& G,
                       const NumericMatrix& g,
                       double scalingFactor) {
  
  // Convert inputs to Armadillo
  mat M_mat = as<mat>(M);        // n x p
  mat G_mat = as<mat>(G);        // n x n
  mat g_mat = as<mat>(g);        // n x k
  int n = M_mat.n_rows;
  int p = M_mat.n_cols;
  int k = g_mat.n_cols;
  
  // Check/regularize G, then invert
  arma::mat L;
  bool spd = arma::chol(L, G_mat);
  if (!spd) {
    Rcpp::Rcout << "G.mat is not psd, calculating nearest psd.\n";
    arma::vec eval;
    arma::mat evec;
    arma::eig_sym(eval, evec,G_mat);
    
    eval.transform([](double xx){ return (xx < 0.0) ? 0.0 : xx; });
    
    G_mat = evec * arma::diagmat(eval) * evec.t();
  }
  
  
  mat Gi;
  bool ok = inv_sympd(Gi, G_mat);
  if (!ok) Gi = inv(G_mat);
  
  arma::mat tmp = Gi * g_mat;                 // n x k
  arma::mat mu  = (scalingFactor * M_mat.t()) * tmp;  // p x k
  
  return List::create(
    Named("mu_matrix")      = mu                              // p x k
   );
}
