// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
DataFrame cpp_calcOHV(const NumericMatrix& Crosses, const List& HBlocks, const NumericMatrix& M, const NumericVector& mu_vec, int nThreads=4) {
  omp_set_num_threads(nThreads); //Sets number of threads for OpenMP
  int numCrosses = Crosses.nrow();
  
  // Result column
  NumericVector OHV_col(numCrosses);
  
  // Convert inputs to Armadillo structures
  rowvec mu = as<rowvec>(mu_vec);
  mat M_mat = as<mat>(M);

  // Calculate total number of columns required for LocalGEBV
  int totalColumns = HBlocks.size();
  mat LocalGEBV(M_mat.n_rows, totalColumns, fill::zeros);

  // Fill LocalGEBV with contributions from each block
  int colOffset = 0;
  for (int i = 0; i < HBlocks.size(); i++) {
    IntegerVector block = HBlocks[i];
    uvec indices = as<uvec>(block) - 1; // Convert to zero-based indexing

    // Boundary checks
    if (indices.max() >= M_mat.n_cols || indices.min() < 0) {
      stop("Index in HBlocks out of range.");
    }

    // Compute block-specific GEBVs
    mat Z_block = M_mat.cols(indices);
    vec mu_block = mu(indices);
    mat blockGEBV = Z_block * mu_block;

    // Assign to LocalGEBV
    LocalGEBV.col(colOffset) = blockGEBV; // Store each blockGEBV as a column
    colOffset += 1;
  }

  // Parallel computation of OHV
#pragma omp parallel for
  for (int x = 0; x < numCrosses; ++x) {
    int P1_index = Crosses(x, 0) - 1;
    int P2_index = Crosses(x, 1) - 1;

    // Validate parent indices
    if (P1_index < 0 || P1_index >= static_cast<int>(M_mat.n_rows) || 
        P2_index < 0 || P2_index >= static_cast<int>(M_mat.n_rows)) {
      continue;
    }

    // Compute OHV for this cross
    rowvec max_values = arma::max(LocalGEBV.row(P1_index), LocalGEBV.row(P2_index));
    OHV_col[x] = arma::sum(max_values);
  }

  return DataFrame::create(
    Named("OHV") = OHV_col
  );
}
