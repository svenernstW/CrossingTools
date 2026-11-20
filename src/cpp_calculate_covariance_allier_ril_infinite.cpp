// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <atomic>
#include <algorithm>
#include <limits>
#include <utility>
#include "parallel_backend.h"   // OpenMP when available, else RcppParallel

using namespace Rcpp;
using namespace arma;

inline double ril_infinite_kernel(double x, double y, int t) {
  double diff = std::abs(x - y);
  double rcf  = 0.5 * (1.0 - std::exp(-2.0 * diff));

  double randomMatingKernel = std::pow(1.0 - rcf, t);
  double infiniteRilKernel = (1.0 - 2.0 * rcf) / (1.0 + 2.0 * rcf);

  return randomMatingKernel * infiniteRilKernel;
}

// [[Rcpp::export]]
SEXP cpp_calculate_covariance_allier_ril_infinite(const NumericMatrix& Crosses,
                                                  const List& genMap,
                                                  const NumericMatrix& M,
                                                  const NumericMatrix& U,
                                                  int t,
                                                  double intensity,
                                                  const NumericVector& weights,
                                                  const NumericVector& p,
                                                  bool covariance = false,
                                                  bool calcindex = false,
                                                  int nThreads = 4) {
  // Portable thread setup
  ct_set_threads(nThreads);
  std::atomic<bool> any_psd(false);

  const arma::uword numCrosses   = Crosses.nrow();
  const arma::uword numMarkers   = M.ncol();
  const arma::uword numTrait     = U.ncol();
  const arma::uword numTraitComb = numTrait * (numTrait + 1) / 2;

  if (t < 0) {
    Rcpp::stop("t must be non-negative");
  }
  if (Crosses.ncol() != 4 || U.nrow() != static_cast<int>(numMarkers) ||
      p.size() != static_cast<R_xlen_t>(numMarkers) ||
      weights.size() != static_cast<R_xlen_t>(numTrait) || numTrait == 0) {
    Rcpp::stop("Invalid Crosses, M, U, p or weights dimensions");
  }

  // Convert inputs to Armadillo
  arma::mat M_mat = as<arma::mat>(M);  // parental genotypes coded 0/2
  arma::mat U_mat = as<arma::mat>(U);  // marker effects (markers x traits)
  const arma::uword nInd = M_mat.n_rows;
  arma::vec weights_vec = as<arma::vec>(weights);
  arma::rowvec p_vec = as<arma::rowvec>(p);

  // Genomic breeding values, without copying the full centred genotype matrix
  arma::mat GEBV = M_mat * U_mat;
  arma::rowvec meanAdjustment = 2.0 * (p_vec * U_mat);
  GEBV.each_row() -= meanAdjustment;

  // Contiguous trait effects for one marker, for the innermost trait loops
  arma::mat U_trans = U_mat.t();

  const arma::uword OFF_EG  = 0;
  const arma::uword OFF_VAR = numTrait;
  const arma::uword OFF_SPV = 2 * numTrait;


  // Read marker positions and identify chromosome boundaries.
  // Only pairs of markers on the same chromosome contribute to the variance.
  std::vector<double> markerPos(numMarkers);
  std::vector<std::pair<arma::uword, arma::uword>> chromosomeRanges;
  arma::uword startIdx = 0;

  for (R_xlen_t c = 0; c < genMap.size(); ++c) {
    NumericMatrix gm = as<NumericMatrix>(genMap[c]);
    const arma::uword n = gm.nrow();
    if (n == 0) continue;
    if (gm.ncol() < 1 || n > numMarkers - startIdx) {
      Rcpp::stop("Invalid genetic map dimensions");
    }

    for (arma::uword j = 0; j < n; ++j) {
      const double pos = gm(j, 0);
      if (!std::isfinite(pos) || (j > 0 && pos < gm(j - 1, 0))) {
        Rcpp::stop("Marker positions must be finite and ordered within chromosomes");
      }
      markerPos[startIdx + j] = pos;
    }

    chromosomeRanges.emplace_back(startIdx, startIdx + n);
    startIdx += n;
  }

  if (startIdx != numMarkers) {
    Rcpp::stop("genMap and M contain different marker counts");
  }

  // Copy parent indices out of Crosses (1-based in R -> 0-based here).
  // Validate before entering the parallel region.
  arma::Mat<int> P(numCrosses, 4);
  for (arma::uword x = 0; x < numCrosses; ++x) {
    for (arma::uword a = 0; a < 4; ++a) {
      const double parent = Crosses(x, a);
      if (!std::isfinite(parent) || parent < 1.0 ||
          parent > static_cast<double>(nInd) ||
          parent != std::floor(parent) ||
          parent > static_cast<double>(std::numeric_limits<int>::max())) {
        Rcpp::stop("Crosses must contain valid integer parent indices");
      }
      P(x, a) = static_cast<int>(parent) - 1;
    }
  }

  // Results (same format as original four-way code)
  arma::mat results1(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results2(numCrosses, numTrait * 3 + 3, arma::fill::zeros);

  auto tri_u_idx_incl = [numTrait](arma::uword ti, arma::uword tj) -> arma::uword {
    return ti * numTrait - (ti * (ti - 1)) / 2 + (tj - ti);
  };

  // Parallel over four-way crosses
  ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
    const arma::uword x = static_cast<arma::uword>(xi);
    const int P1 = P(x, 0), P2 = P(x, 1), P3 = P(x, 2), P4 = P(x, 3);

    // Expected genetic values (mean of the four parents)
    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      const double eG = 0.25 * (GEBV(P1, ti) + GEBV(P2, ti) +
                                GEBV(P3, ti) + GEBV(P4, ti));
      results2(x, OFF_EG + ti) = eG;
      results2(x, OFF_SPV + ti) = eG;
    }

    // Accumulate contributions from segregating markers on each chromosome.
    for (const auto& chr : chromosomeRanges) {
      const arma::uword startC = chr.first;
      const arma::uword endC = chr.second;  // exclusive end

      std::vector<arma::uword> differing;
      std::vector<double> a12, a34, between;
      differing.reserve(endC - startC);
      a12.reserve(endC - startC);
      a34.reserve(endC - startC);
      between.reserve(endC - startC);

      // Exclude markers for which all four parents have identical genotypes.
      for (arma::uword mi = startC; mi < endC; ++mi) {
        const double m1 = M_mat(P1, mi);
        const double m2 = M_mat(P2, mi);
        const double m3 = M_mat(P3, mi);
        const double m4 = M_mat(P4, mi);
        if (m1 == m2 && m1 == m3 && m1 == m4) continue;

        differing.push_back(mi);
        a12.push_back(m1 - m2);
        a34.push_back(m3 - m4);
        between.push_back(m1 + m2 - m3 - m4);
      }

      const std::size_t mlen = differing.size();
      if (mlen == 0) continue;

      for (std::size_t a = 0; a < mlen; ++a) {
        const arma::uword gi = differing[a];

        for (std::size_t b = a; b < mlen; ++b) {
          const arma::uword gj = differing[b];

          // Two within-block contrasts and one across-block contrast.
          const double phi =
            2.0 * a12[a] * a12[b]
          + 2.0 * a34[a] * a34[b]
          + between[a] * between[b];
          if (phi == 0.0) continue;

          // Calculate the exact infinite-RIL kernel directly for this marker pair.
          const double kernel = ril_infinite_kernel(markerPos[gi], markerPos[gj], t);

          // Both the 1/4 outer factor and the 1/4 in cB are included here.
          const double Dcomb = 0.0625 * kernel * phi;
          if (Dcomb == 0.0) continue;

          if (!covariance) {
            // Trait variances only; no off-diagonal trait computations.
            const double mult = (a == b) ? Dcomb : 2.0 * Dcomb;
            for (arma::uword ti = 0; ti < numTrait; ++ti) {
              const arma::uword kcol = tri_u_idx_incl(ti, ti);
              results1(x, kcol) += mult * U_trans(ti, gi) * U_trans(ti, gj);
            }
          } else {
            // Only the upper triangle of the trait covariance matrix.
            for (arma::uword ti = 0; ti < numTrait; ++ti) {
              const double Ui_ti = U_trans(ti, gi);
              const double Uj_ti = U_trans(ti, gj);
              for (arma::uword tj = ti; tj < numTrait; ++tj) {
                const arma::uword kcol = tri_u_idx_incl(ti, tj);
                if (a == b) {
                  results1(x, kcol) += Dcomb * Ui_ti * U_trans(tj, gi);
                } else {
                  results1(x, kcol) += Dcomb * (
                    Ui_ti * U_trans(tj, gj) +
                      Uj_ti * U_trans(tj, gi));
                }
              }
            }
          }
        }
      }
    }

    // Expected value, segregation variance, and superior progeny value
    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      const arma::uword kdiag = tri_u_idx_incl(ti, ti);
      const double var = results1(x, kdiag);
      const double eG = results2(x, OFF_EG + ti);
      results2(x, OFF_VAR + ti) = var;
      results2(x, OFF_SPV + ti) =
        eG + intensity * std::sqrt(std::max(0.0, var));
    }
  });

  // If requested, convert upper-triangle values to symmetric trait matrices.
  if (covariance) {
    std::vector<arma::mat> covs(numCrosses);

    ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
      const arma::uword x = static_cast<arma::uword>(xi);
      arma::mat V(numTrait, numTrait, arma::fill::zeros);
      for (arma::uword ti = 0; ti < numTrait; ++ti) {
        for (arma::uword tj = ti; tj < numTrait; ++tj) {
          const arma::uword k = tri_u_idx_incl(ti, tj);
          const double value = results1(x, k);
          V(ti, tj) = value;
          V(tj, ti) = value;
        }
      }
      covs[x] = V;

      if (calcindex) {
        arma::mat L;
        if (!arma::chol(L, V)) {
          any_psd.store(true, std::memory_order_relaxed);
          arma::vec eval;
          arma::mat evec;
          arma::eig_sym(eval, evec, V);
          eval.transform([](double val) { return (val < 0.0) ? 0.0 : val; });
          V = evec * arma::diagmat(eval) * evec.t();
        }
        const double sigma_gains = arma::as_scalar(weights_vec.t() * V * weights_vec);
        const double index = arma::as_scalar(
          results2.row(x).cols(0, numTrait - 1) * weights_vec);
        const double spvi = index + intensity * std::sqrt(std::max(0.0, sigma_gains));
        results2(x, 3 * numTrait)     = index;
        results2(x, 3 * numTrait + 1) = sigma_gains;
        results2(x, 3 * numTrait + 2) = spvi;
      }
    });

    Rcpp::List out(numCrosses);
    for (arma::uword x = 0; x < numCrosses; ++x) out[x] = Rcpp::wrap(covs[x]);
    return Rcpp::List::create(
      Rcpp::Named("cross_values") = results2,
      Rcpp::Named("covariances")  = out,
      Rcpp::Named("check_psd")   = any_psd.load(std::memory_order_relaxed)
    );
  }

  return Rcpp::wrap(results2);
}
