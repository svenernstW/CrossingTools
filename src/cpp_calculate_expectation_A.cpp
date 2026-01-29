// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <atomic>

#include "parallel_backend.h"  // OpenMP when available, else RcppParallel

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
SEXP cpp_calculate_expectation_A(const NumericMatrix& Crosses,
                                  const NumericMatrix& M,
                                  const NumericMatrix& U,
                                  const NumericVector& weights,
                                  bool calcindex = false,
                                  int nThreads = 4) {
  // portable thread setup
  ct_set_threads(nThreads);
  
  const arma::uword numCrosses   = Crosses.nrow();
  const arma::uword numTrait     = U.ncol();
  
  // Convert inputs to Armadillo
  arma::mat M_mat = as<arma::mat>(M);
  arma::mat U_mat    = as<arma::mat>(U);      // (numMarkers × numTrait)
  arma::vec weights_vec = as<arma::vec>(weights); // length == numTrait
  
  const arma::uword OFF_EG    = 0;
  
  // Copy parent indices out of Crosses (1-based in R → 0-based here)
  arma::Col<int> P1_idx(numCrosses), P2_idx(numCrosses);
  for (arma::uword x = 0; x < numCrosses; ++x) {
    P1_idx[x] = static_cast<int>(Crosses(x, 0)) - 1;
    P2_idx[x] = static_cast<int>(Crosses(x, 1)) - 1;
  }
  
  // Results
  arma::mat results2(numCrosses, numTrait+1, arma::fill::zeros);
  
  // Parallel over crosses
  ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
    R_xlen_t x = static_cast<R_xlen_t>(xi);
    const int P1 = P1_idx[x];
    const int P2 = P2_idx[x];
    const arma::uword nInd = M_mat.n_rows;
    if (P1 < 0 || P2 < 0 || P1 >= (int)nInd || P2 >= (int)nInd) return;
    
    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      
      double G1 = arma::dot(M_mat.row(P1), U_mat.col(ti));
      double G2 = arma::dot(M_mat.row(P2), U_mat.col(ti));
      double eG = 0.5 * (G1 + G2);
      
      results2(x, OFF_EG   + ti) = eG;
    }
    
  });
  
  if (calcindex) {
    
    // Offsets for the 6 scalar outputs (indices/variances/SPV for A-only and A+D)
    const arma::uword OFF_IDX_A     = 1 * numTrait + 0;
    
    
    ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
      arma::uword x = static_cast<arma::uword>(xi);
      
      
      
      double index_A = arma::as_scalar(results2.row(x).cols(0, numTrait - 1) * weights_vec);
      
      results2(x, OFF_IDX_A)    = index_A;
      
        
      
    });
  }
  
  
  
  return Rcpp::wrap(results2);
}
