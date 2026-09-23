// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <atomic>

#include "parallel_backend.h"  // OpenMP when available, else RcppParallel

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
SEXP cpp_calculate_expectation_AD(const NumericMatrix& Crosses,
                                    const NumericMatrix& Hap1,
                                    const NumericMatrix& Hap2,
                                    const NumericMatrix& U,
                                    const NumericMatrix& D,
                                    const NumericVector& weights,
                                    const NumericVector& p,
                                    bool calcindex = false,
                                    int nThreads = 4) {
  // portable thread setup
  ct_set_threads(nThreads);

  const arma::uword numCrosses   = Crosses.nrow();
  const arma::uword numTrait     = U.ncol();

  // Convert inputs to Armadillo
  arma::mat Hap1_mat = as<arma::mat>(Hap1);   // (n_individuals × numMarkers)
  arma::mat Hap2_mat = as<arma::mat>(Hap2);   // (n_individuals × numMarkers)
  arma::mat M_mat    = Hap1_mat + Hap2_mat;   // dosage 0..2
  arma::mat U_mat    = as<arma::mat>(U);      // average substitution effects alpha
  arma::mat D_mat    = as<arma::mat>(D);      // dominance effects d
  arma::vec weights_vec = as<arma::vec>(weights);

  // Allele frequencies in the reference population
  arma::vec p_vec = as<arma::vec>(p);

  // Centred additive genotype code W = M - 2p
  arma::mat W_mat = M_mat;
  W_mat.each_row() -= 2.0 * p_vec;

  // Genomic breeding values
  arma::mat GEBV =
    W_mat * U_mat;
  const arma::uword OFF_EG    = 0;
  const arma::uword OFF_ETG   = numTrait;

  // Copy parent indices out of Crosses (1-based in R → 0-based here)
  arma::Col<int> P1_idx(numCrosses), P2_idx(numCrosses);
  for (arma::uword x = 0; x < numCrosses; ++x) {
    P1_idx[x] = static_cast<int>(Crosses(x, 0)) - 1;
    P2_idx[x] = static_cast<int>(Crosses(x, 1)) - 1;
  }

  // Results
  arma::mat results2(numCrosses, numTrait * 2 + 2, arma::fill::zeros);

  // Parallel over crosses
  ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
    R_xlen_t x = static_cast<R_xlen_t>(xi);
    const int P1 = P1_idx[x];
    const int P2 = P2_idx[x];
    const arma::uword nInd = M_mat.n_rows;
    if (P1 < 0 || P2 < 0 || P1 >= (int)nInd || P2 >= (int)nInd) return;

    const arma::rowvec W1 =
      W_mat.row(P1);

    const arma::rowvec W2 =
      W_mat.row(P2);

    // Expected statistical dominance-deviation code
    const arma::colvec meanD =
      (-0.5 * (W1 % W2)).t();

    for (arma::uword ti = 0; ti < numTrait; ++ti) {

      const double eG =
        0.5 * (
            GEBV(P1, ti) +
              GEBV(P2, ti)
        );

      const double eD =
        arma::dot(
          D_mat.col(ti),
          meanD
        );

      const double eTG =
        eG + eD;

      results2(x, OFF_EG  + ti) =
        eG;

      results2(x, OFF_ETG + ti) =
        eTG;
    }

  });

  if (calcindex) {

    const arma::uword OFF_IDX_A =
      2 * numTrait + 0;

    const arma::uword OFF_IDX_T =
      2 * numTrait + 1;

    ct_parallel_for(
      0,
      static_cast<int>(numCrosses),
      [&](int xi) {

        const arma::uword x =
          static_cast<arma::uword>(xi);

        // Expected breeding-value index
        const double index_A =
          arma::as_scalar(
            results2.row(x)
                    .cols(OFF_EG, OFF_EG + numTrait - 1) *
            weights_vec
          );

        // Expected total-genetic-value index
        const double index_T =
          arma::as_scalar(
            results2.row(x)
                    .cols(OFF_ETG, OFF_ETG + numTrait - 1) *
            weights_vec
          );

        results2(x, OFF_IDX_A) =
          index_A;

        results2(x, OFF_IDX_T) =
          index_T;
      }
    );
  }


  return Rcpp::wrap(results2);
}
