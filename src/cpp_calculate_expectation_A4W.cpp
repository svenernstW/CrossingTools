// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <atomic>

#include "parallel_backend.h"  // OpenMP when available, else RcppParallel

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
SEXP cpp_calculate_expectation_A4W(const NumericMatrix& Crosses,
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
  //Precompute GEBV
  arma::mat GEBV = M_mat * U_mat;  // (nInd × numTrait)

  const arma::uword OFF_EG    = 0;

  // Copy parent indices out of Crosses (1-based in R → 0-based here)
  arma::Col<int> P1_idx(numCrosses), P2_idx(numCrosses), P3_idx(numCrosses), P4_idx(numCrosses);
  for (arma::uword x = 0; x < numCrosses; ++x) {
    P1_idx[x] = static_cast<int>(Crosses(x, 0)) - 1;
    P2_idx[x] = static_cast<int>(Crosses(x, 1)) - 1;
    P3_idx[x] = static_cast<int>(Crosses(x, 2)) - 1;
    P4_idx[x] = static_cast<int>(Crosses(x, 3)) - 1;
  }

  // Results
  arma::mat results2(numCrosses, numTrait+1, arma::fill::zeros);

  // Parallel over crosses
  ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
    R_xlen_t x = static_cast<R_xlen_t>(xi);
    const int P1 = P1_idx[x];
    const int P2 = P2_idx[x];
    const int P3 = P3_idx[x];
    const int P4 = P4_idx[x];
    const arma::uword nInd = M_mat.n_rows;
    if (P1 < 0 || P2 < 0 || P1 >= (int)nInd || P2 >= (int)nInd || P3 < 0 || P4 < 0 || P3 >= (int)nInd || P4 >= (int)nInd) return;

    for (arma::uword ti = 0; ti < numTrait; ++ti) {

      double G1 = GEBV(P1, ti);
      double G2 = GEBV(P2, ti);
      double G3 = GEBV(P3, ti);
      double G4 = GEBV(P4, ti);
      double eG = 0.25 * (G1 + G2 + G3 + G4);

      results2(x, OFF_EG   + ti) = eG;
    }

  });

  if (calcindex) {
    //Precompute GEBV
    arma::vec GEBVindex = GEBV * weights_vec;   // (nInd)

    // Offsets for the 6 scalar outputs (indices/variances/SPV for A-only and A+D)
    const arma::uword OFF_IDX_A     = 1 * numTrait + 0;


    ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
      arma::uword x = static_cast<arma::uword>(xi);
      int P1 = P1_idx[x];
      int P2 = P2_idx[x];
      int P3 = P3_idx[x];
      int P4 = P4_idx[x];
      arma::uword nInd = M_mat.n_rows;
      if (P1 < 0 || P2 < 0 || P1 >= (int)nInd || P2 >= (int)nInd || P3 < 0 || P4 < 0 || P3 >= (int)nInd || P4 >= (int)nInd) return;

      double eGindex = 0.25 * (GEBVindex(P1) + GEBVindex(P2)+GEBVindex(P3) + GEBVindex(P4));
      results2(x, OFF_IDX_A) = eGindex;



    });
  }



  return Rcpp::wrap(results2);
}
