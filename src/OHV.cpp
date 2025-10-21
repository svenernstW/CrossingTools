// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include "parallel_backend.h"  // OpenMP when available, else RcppParallel

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
DataFrame cpp_calcOHV(const NumericMatrix& Crosses,
                      const List& HBlocks,
                      const NumericMatrix& M,
                      const NumericVector& mu_vec,
                      int nThreads = 4) {
  // portable thread setup
  ct_set_threads(nThreads);

  const int numCrosses = Crosses.nrow();

  // Result column
  NumericVector OHV_col(numCrosses);

  // Convert inputs to Armadillo
  rowvec mu_row = as<rowvec>(mu_vec);  // length = n_mark
  mat    M_mat  = as<mat>(M);          // (n_ind x n_mark)

  // Precompute per-block local GEBV columns
  const int nInd = static_cast<int>(M_mat.n_rows);
  const int nMark = static_cast<int>(M_mat.n_cols);
  const int nBlocks = HBlocks.size();

  mat LocalGEBV(M_mat.n_rows, nBlocks, fill::zeros);

  int colOffset = 0;
  for (int i = 0; i < nBlocks; ++i) {
    // 1-based indices from R
    IntegerVector block = HBlocks[i];
    if (block.size() == 0) {
      stop("HBlocks[[%d]] is empty.", i + 1);
    }
    // Validate range (1..nMark)
    int bmin = block[0], bmax = block[0];
    for (int j = 0; j < block.size(); ++j) {
      int v = block[j];
      if (v < 1 || v > nMark) {
        stop("Index in HBlocks[[%d]] out of range [1, %d].", i + 1, nMark);
      }
      if (v < bmin) bmin = v;
      if (v > bmax) bmax = v;
    }
    // Convert to 0-based uvec
    uvec indices(block.size());
    for (int j = 0; j < block.size(); ++j) indices[j] = static_cast<uword>(block[j] - 1);

    // Compute block-specific GEBV = Z_block * mu_block
    mat Z_block = M_mat.cols(indices);                                 // (n_ind x |block|)
    vec mu_block = conv_to<vec>::from(mu_row.elem(indices));           // (|block| x 1)
    vec blockGEBV = Z_block * mu_block;                                // (n_ind x 1)

    LocalGEBV.col(colOffset) = blockGEBV; // store each block as a column
    colOffset += 1;
  }

  // Parallel over crosses
  ct_parallel_for(0, numCrosses, [&](int x) {
    int P1_index = static_cast<int>(Crosses(x, 0)) - 1;
    int P2_index = static_cast<int>(Crosses(x, 1)) - 1;

    // Validate parent indices
    if (P1_index < 0 || P1_index >= nInd || P2_index < 0 || P2_index >= nInd) {
      // leave OHV_col[x] as 0
      return;
    }

    // OHV: sum over blocks of max(LocalGEBV[P1, b], LocalGEBV[P2, b])
    rowvec max_values = arma::max(LocalGEBV.row(P1_index), LocalGEBV.row(P2_index));
    OHV_col[x] = arma::sum(max_values);
  });

  return DataFrame::create(Named("OHV") = OHV_col);
}
