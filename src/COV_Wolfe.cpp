// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include "parallel_backend.h"  // OpenMP when available, else RcppParallel

using namespace Rcpp;
using namespace arma;

// Inline function to calculate recombination fraction
inline double qjk(double x, double y) {
  double diff = std::abs(x - y);
  double rcf  = 0.5 * (1 - std::exp(-2 * diff));
  return 1 - 2 * rcf;
}

// [[Rcpp::export]]
SEXP cpp_calculate_covariance_wolfe(const NumericMatrix& Crosses,
                                    const List& genMap,
                                    const NumericMatrix& Hap1,
                                    const NumericMatrix& Hap2,
                                    const NumericMatrix& U,
                                    const NumericMatrix& D,
                                    double intensity,
                                    const NumericVector& gains,
                                    bool covariance = false,
                                    bool calcgains = false,
                                    int nThreads = 4) {
  // portable thread setup
  ct_set_threads(nThreads);

  const arma::uword numCrosses   = Crosses.nrow();
  const arma::uword numMarkers   = Hap1.ncol();
  const arma::uword numTrait     = U.ncol();
  const arma::uword numTraitComb = numTrait * (numTrait + 1) / 2;

  // Convert inputs to Armadillo
  arma::mat Hap1_mat = as<arma::mat>(Hap1);   // (n_individuals × numMarkers)
  arma::mat Hap2_mat = as<arma::mat>(Hap2);   // (n_individuals × numMarkers)
  arma::mat M_mat    = Hap1_mat + Hap2_mat;   // dosage 0..2
  arma::mat U_mat    = as<arma::mat>(U);      // (numMarkers × numTrait)
  arma::mat D_mat    = as<arma::mat>(D);      // (numMarkers × numTrait)
  arma::vec gains_vec = as<arma::vec>(gains); // length == numTrait

  const arma::uword OFF_EG    = 0;
  const arma::uword OFF_ETG   = numTrait;
  const arma::uword OFF_VAR_A = 2 * numTrait;
  const arma::uword OFF_SPV   = 3 * numTrait;
  const arma::uword OFF_VAR_D = 4 * numTrait;
  const arma::uword OFF_TSPV  = 5 * numTrait;

  // Precompute chromosome ranges
  std::vector<std::pair<int, int>> chromosomeRanges;
  {
    int startIdx = 0;
    for (int i = 0; i < genMap.size(); ++i) {
      NumericMatrix genMapMatrix = as<NumericMatrix>(genMap[i]);
      int n = genMapMatrix.nrow();
      chromosomeRanges.emplace_back(startIdx, startIdx + n - 1);
      startIdx += n;
    }
  }

  // Precompute recombination correlation RC(i,j) = 1 - 2r_ij
  arma::mat RC(numMarkers, numMarkers, fill::zeros);
  {
    int startIdx = 0;
    for (int i = 0; i < genMap.size(); ++i) {
      NumericMatrix genMapMatrix = as<NumericMatrix>(genMap[i]);
      int n = genMapMatrix.nrow();
      for (int j = 0; j < n; ++j) {
        for (int k = j; k < n; ++k) {
          double value = qjk(genMapMatrix(j, 0), genMapMatrix(k, 0));
          RC(startIdx + j, startIdx + k) = value;
          RC(startIdx + k, startIdx + j) = value;
        }
      }
      startIdx += n;
    }
  }

  // Copy parent indices out of Crosses (1-based in R → 0-based here)
  arma::Col<int> P1_idx(numCrosses), P2_idx(numCrosses);
  for (arma::uword x = 0; x < numCrosses; ++x) {
    P1_idx[x] = static_cast<int>(Crosses(x, 0)) - 1;
    P2_idx[x] = static_cast<int>(Crosses(x, 1)) - 1;
  }

  // Results
  arma::mat results1A(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results1D(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results2(numCrosses, numTrait * 6 + 6, arma::fill::zeros);

  // Helper to map (ti, tj) with 0 ≤ ti ≤ tj < numTrait to column index
  auto tri_u_idx_incl = [numTrait](arma::uword ti, arma::uword tj) -> arma::uword {
    return ti * numTrait - (ti * (ti - 1)) / 2 + (tj - ti);
  };

  // Parallel over crosses
  ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
    R_xlen_t x = static_cast<R_xlen_t>(xi);
    const int P1 = P1_idx[x];
    const int P2 = P2_idx[x];
    const arma::uword nInd = M_mat.n_rows;
    if (P1 < 0 || P2 < 0 || P1 >= (int)nInd || P2 >= (int)nInd) return;

    // Markers that are heterozygous in exactly one parent (dosage == 1 in either parent)
    arma::uvec differingMarkers = arma::find(
      (M_mat.row(P1) == 1.0) + (M_mat.row(P2) == 1.0)
    );

    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      for (arma::uword tj = ti; tj < numTrait; ++tj) {
        if (!covariance && ti != tj) continue;

        double SigmaSqAP1P2 = 0.0;
        double SigmaSqDP1P2 = 0.0;

        for (const auto& chrRange : chromosomeRanges) {
          int startC = chrRange.first;
          int endC   = chrRange.second;

          arma::uvec chrDiff = differingMarkers( arma::find(
            (differingMarkers >= (arma::uword)startC) %
              (differingMarkers <= (arma::uword)endC)
          ));
          if (chrDiff.n_elem == 0) continue;

          for (uword i = 0; i < chrDiff.n_elem; ++i) {
            int marker_i = (int)chrDiff[i];
            for (uword j = i; j < chrDiff.n_elem; ++j) {
              int marker_j = (int)chrDiff[j];

              // Recombination correlation
              double rcf_ij = RC(marker_i, marker_j);

              // Build small 2×2 haplotype matrices per parent (markers i,j)
              arma::uword row_idx_P1 = static_cast<arma::uword>(P1);
              arma::uword row_idx_P2 = static_cast<arma::uword>(P2);
              arma::uvec cols(2); cols(0) = (arma::uword)marker_i; cols(1) = (arma::uword)marker_j;

              arma::rowvec Hap1_P1 = Hap1_mat(arma::uvec{row_idx_P1}, cols);
              arma::rowvec Hap2_P1 = Hap2_mat(arma::uvec{row_idx_P1}, cols);
              arma::rowvec Hap1_P2 = Hap1_mat(arma::uvec{row_idx_P2}, cols);
              arma::rowvec Hap2_P2 = Hap2_mat(arma::uvec{row_idx_P2}, cols);

              arma::rowvec P1_means = 0.5 * (Hap1_P1 + Hap2_P1);
              arma::rowvec P2_means = 0.5 * (Hap1_P2 + Hap2_P2);

              arma::mat Hap_P1 = join_cols(Hap1_P1, Hap2_P1);
              arma::mat Hap_P2 = join_cols(Hap1_P2, Hap2_P2);

              arma::mat D1 = 0.5 * Hap_P1.t() * Hap_P1 - P1_means.t() * P1_means;
              arma::mat D2 = 0.5 * Hap_P2.t() * Hap_P2 - P2_means.t() * P2_means;

              // DGen for the pair (i,j)
              double DGen_ij = (rcf_ij * D1.at(0,1)) + (rcf_ij * D2.at(0,1));

              double contribA = U_mat(marker_j, tj) * DGen_ij * U_mat(marker_i, ti);
              SigmaSqAP1P2 += (i == j) ? contribA : 2 * contribA;

              double contribD = D_mat(marker_j, tj) * (DGen_ij * DGen_ij) * D_mat(marker_i, ti);
              SigmaSqDP1P2 += (i == j) ? contribD : 2 * contribD;
            }
          }
        }

        const arma::uword k = tri_u_idx_incl(ti, tj);
        results1A(x, k) = SigmaSqAP1P2;
        results1D(x, k) = SigmaSqDP1P2;

        if (ti == tj) {
          double G1 = arma::dot(M_mat.row(P1), U_mat.col(ti));
          double G2 = arma::dot(M_mat.row(P2), U_mat.col(ti));
          double eG = 0.5 * (G1 + G2);

          arma::rowvec pk = M_mat.row(P1);
          arma::rowvec qk = 2.0 - pk;
          arma::rowvec ykl = pk - M_mat.row(P2);
          const arma::colvec v1 = (pk - qk - ykl).t();
          const arma::colvec v2 = (2.0 * (pk % qk) + (ykl % (pk - qk))).t();

          double eTG = arma::dot(U_mat.col(ti), v1) + arma::dot(D_mat.col(ti), v2);

          results2(x, OFF_EG   + ti) = eG;
          results2(x, OFF_VAR_A+ ti) = SigmaSqAP1P2;
          results2(x, OFF_SPV  + ti) = eG + intensity * std::sqrt(SigmaSqAP1P2);
          results2(x, OFF_VAR_D+ ti) = SigmaSqDP1P2;
          results2(x, OFF_ETG  + ti) = eTG;
          results2(x, OFF_TSPV + ti) = eTG + (std::sqrt(SigmaSqAP1P2) + std::sqrt(SigmaSqDP1P2));
        }
      }
    }
  });

  if (covariance) {
    std::vector<arma::mat> covsA(numCrosses);
    std::vector<arma::mat> covsD(numCrosses);

    // Offsets for the 6 scalar outputs (indices/variances/SPV for A-only and A+D)
    const arma::uword OFF_IDX_A     = 6 * numTrait + 0;
    const arma::uword OFF_VARIDX_A  = 6 * numTrait + 1;
    const arma::uword OFF_SPVI_A    = 6 * numTrait + 2;
    const arma::uword OFF_IDX_AD    = 6 * numTrait + 3;
    const arma::uword OFF_VARIDX_AD = 6 * numTrait + 4;
    const arma::uword OFF_SPVI_AD   = 6 * numTrait + 5;

    ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
      arma::uword x = static_cast<arma::uword>(xi);

      // Build symmetric covariance matrices GA (additive) and GD (dominance)
      arma::mat GA(numTrait, numTrait, arma::fill::zeros);
      arma::mat GD(numTrait, numTrait, arma::fill::zeros);
      for (arma::uword ti = 0; ti < numTrait; ++ti) {
        for (arma::uword tj = ti; tj < numTrait; ++tj) {
          const arma::uword k =
            (ti * numTrait) - (ti * (ti - 1)) / 2 + (tj - ti);
          const double va = results1A(x, k);
          const double vd = results1D(x, k);
          GA(ti, tj) = GA(tj, ti) = va;
          GD(ti, tj) = GD(tj, ti) = vd;
        }
      }
      covsA[x] = GA;
      covsD[x] = GD;

      if (calcgains) {
        arma::mat VA = GA;
        arma::mat L;
        if (!arma::chol(L, VA)) {
          arma::vec eval; arma::mat evec;
          arma::eig_sym(eval, evec, VA);
          double tol = std::max(1e-12, 1e-8 * eval.max());
          for (double &lam : eval) if (lam < tol) lam = tol;
          VA = evec * arma::diagmat(eval) * evec.t();
          arma::chol(L, VA);
        }
        arma::vec wA = arma::solve(arma::trimatu(L.t()),
                                   arma::solve(arma::trimatl(L), gains_vec));

        double var_index_A = arma::as_scalar(wA.t() * GA * wA);
        var_index_A = std::max(0.0, var_index_A); // numeric safety
        double index_A = arma::as_scalar(results2.row(x).cols(0, numTrait - 1) * wA);
        double spvi_A  = index_A + intensity * std::sqrt(var_index_A);

        results2(x, OFF_IDX_A)    = index_A;
        results2(x, OFF_VARIDX_A) = var_index_A;
        results2(x, OFF_SPVI_A)   = spvi_A;

        arma::mat VAD = GA + GD;
        if (!arma::chol(L, VAD)) {
          arma::vec eval; arma::mat evec;
          arma::eig_sym(eval, evec, VAD);
          double tol = std::max(1e-12, 1e-8 * eval.max());
          for (double &lam : eval) if (lam < tol) lam = tol;
          VAD = evec * arma::diagmat(eval) * evec.t();
          arma::chol(L, VAD);
        }
        arma::vec wAD = arma::solve(arma::trimatu(L.t()),
                                    arma::solve(arma::trimatl(L), gains_vec));

        double var_index_AD = arma::as_scalar(wAD.t() * VAD * wAD);
        var_index_AD = std::max(0.0, var_index_AD);
        double index_AD = arma::as_scalar(results2.row(x).cols(numTrait, 2 * numTrait - 1) * wAD);
        double spvi_AD  = index_AD + intensity * std::sqrt(var_index_AD);

        results2(x, OFF_IDX_AD)    = index_AD;
        results2(x, OFF_VARIDX_AD) = var_index_AD;
        results2(x, OFF_SPVI_AD)   = spvi_AD;
      }
    });

    Rcpp::List out_A(numCrosses), out_D(numCrosses);
    for (arma::uword x = 0; x < numCrosses; ++x) {
      out_A[x] = Rcpp::wrap(covsA[x]);
      out_D[x] = Rcpp::wrap(covsD[x]);
    }

    return Rcpp::List::create(
      Rcpp::Named("cross_values") = results2,
      Rcpp::Named("covA")         = out_A,
      Rcpp::Named("covD")         = out_D
    );
  }

  return Rcpp::wrap(results2);
}
