// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <atomic>

#include "parallel_backend.h"  // OpenMP when available, else RcppParallel

using namespace Rcpp;
using namespace arma;

// Inline function to calculate recombination fraction
inline double qjk_wolfe(double x, double y) {
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
                                    const NumericVector& weights,
                                    bool covariance = false,
                                    bool calcindex = false,
                                    int nThreads = 4) {
  // portable thread setup
  ct_set_threads(nThreads);
  std::atomic<bool> any_psd(false);

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
  arma::vec weights_vec = as<arma::vec>(weights); // length == numTrait
  arma::mat GEBV = M_mat * U_mat;  // (nInd × numTrait)

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
          double value = qjk_wolfe(genMapMatrix(j, 0), genMapMatrix(k, 0));
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
    const arma::uword x = static_cast<arma::uword>(xi);
    const int P1 = P1_idx[x];
    const int P2 = P2_idx[x];
    const arma::uword nInd = M_mat.n_rows;
    if (P1 < 0 || P2 < 0 || P1 >= (int)nInd || P2 >= (int)nInd) return;

    // -------------------------
    // (0) Precompute eG + eTG once per trait (diagonal outputs only)
    // -------------------------
    // compute these ONCE per cross
    const arma::rowvec pk  = M_mat.row(P1);
    const arma::rowvec qk  = 2.0 - pk;
    const arma::rowvec ykl = pk - M_mat.row(P2);
    const arma::colvec v1  = (pk - qk - ykl).t();
    const arma::colvec v2  = (2.0 * (pk % qk) + (ykl % (pk - qk))).t();

    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      const double eG  = 0.5 * (GEBV(P1, ti) + GEBV(P2, ti));
      results2(x, OFF_EG + ti) = eG;

      const double eTG = arma::dot(U_mat.col(ti), v1) + arma::dot(D_mat.col(ti), v2);
      results2(x, OFF_ETG + ti) = eTG;
    }


    // -------------------------
    // (1) Build per-chromosome list of "eligible" markers
    // -------------------------
    struct ChrWork {
      std::vector<arma::uword> g; // global marker index
    };
    std::vector<ChrWork> work(chromosomeRanges.size());

    bool any_diff = false;

    for (std::size_t cidx = 0; cidx < chromosomeRanges.size(); ++cidx) {
      const int startC = chromosomeRanges[cidx].first;
      const int endC   = chromosomeRanges[cidx].second;

      auto &w = work[cidx];
      w.g.clear();
      w.g.reserve(256);

      for (int gi = startC; gi <= endC; ++gi) {
        const double m1 = M_mat(P1, (arma::uword)gi);
        const double m2 = M_mat(P2, (arma::uword)gi);

        // Your original criterion was:
        //   (M(P1)==1) OR (M(P2)==1)
        // If you truly mean "exactly one parent het", use XOR:
        //   ((m1==1.0) != (m2==1.0))
        if (!(m1 == 1.0 || m2 == 1.0)) continue;

        any_diff = true;
        w.g.push_back((arma::uword)gi);
      }
    }

    // If nothing contributes, fill diagonal outputs and return
    if (!any_diff) {
      for (arma::uword ti = 0; ti < numTrait; ++ti) {
        const double eG  = results2(x, OFF_EG + ti);
        const double eTG = results2(x, OFF_ETG + ti);

        results2(x, OFF_VAR_A + ti) = 0.0;
        results2(x, OFF_VAR_D + ti) = 0.0;
        results2(x, OFF_SPV   + ti) = eG;
        results2(x, OFF_TSPV  + ti) = eTG;
      }
      return;
    }

    // -------------------------
    // (2) Accumulate variances/covariances:
    //     compute DGen_ij ONCE per marker pair, reuse across traits
    // -------------------------
    if (!covariance) {
      // diagonal only: fill (ti,ti)
      for (std::size_t cidx = 0; cidx < work.size(); ++cidx) {
        const auto &w = work[cidx];
        const std::size_t k = w.g.size();
        if (k == 0) continue;

        for (std::size_t ii = 0; ii < k; ++ii) {
          const arma::uword mi = w.g[ii];

          // parent hap values at marker i
          const double h1i_p1 = Hap1_mat(P1, mi);
          const double h2i_p1 = Hap2_mat(P1, mi);
          const double h1i_p2 = Hap1_mat(P2, mi);
          const double h2i_p2 = Hap2_mat(P2, mi);

          for (std::size_t jj = ii; jj < k; ++jj) {
            const arma::uword mj = w.g[jj];

            const double mult = (ii == jj) ? 1.0 : 2.0;

            // recombination correlation
            const double rcf = RC(mi, mj);

            // parent hap values at marker j
            const double h1j_p1 = Hap1_mat(P1, mj);
            const double h2j_p1 = Hap2_mat(P1, mj);
            const double h1j_p2 = Hap1_mat(P2, mj);
            const double h2j_p2 = Hap2_mat(P2, mj);

            // Off-diagonal element of:
            //   D1 = 0.5 * H^T H - mean mean^T
            // where H has rows hap1 and hap2, and mean = 0.5*(hap1+hap2)
            auto offdiag_D = [](double a_i, double a_j, double b_i, double b_j) -> double {
              const double mean_i = 0.5 * (a_i + b_i);
              const double mean_j = 0.5 * (a_j + b_j);
              return 0.5 * (a_i * a_j + b_i * b_j) - (mean_i * mean_j);
            };

            const double d1 = offdiag_D(h1i_p1, h1j_p1, h2i_p1, h2j_p1);
            const double d2 = offdiag_D(h1i_p2, h1j_p2, h2i_p2, h2j_p2);

            const double DGen = rcf * (d1 + d2);
            const double DGen2 = DGen * DGen;

            // reuse DGen for all traits (diagonal only)
            for (arma::uword ti = 0; ti < numTrait; ++ti) {
              const arma::uword kdiag = tri_u_idx_incl(ti, ti);

              const double contribA = U_mat(mj, ti) * DGen  * U_mat(mi, ti);
              const double contribD = D_mat(mj, ti) * DGen2 * D_mat(mi, ti);

              results1A(x, kdiag) += mult * contribA;
              results1D(x, kdiag) += mult * contribD;
            }
          }
        }
      }

      // write trait-wise outputs (diagonal)
      for (arma::uword ti = 0; ti < numTrait; ++ti) {
        const arma::uword kdiag = tri_u_idx_incl(ti, ti);
        const double varA = results1A(x, kdiag);
        const double varD = results1D(x, kdiag);

        const double eG  = results2(x, OFF_EG  + ti);
        const double eTG = results2(x, OFF_ETG + ti);

        results2(x, OFF_VAR_A + ti) = varA;
        results2(x, OFF_VAR_D + ti) = varD;
        results2(x, OFF_SPV   + ti) = eG  + intensity * std::sqrt(varA);
        results2(x, OFF_TSPV  + ti) = eTG + (std::sqrt(varA) + std::sqrt(varD));
      }

    } else {
      // full covariance: fill all (ti,tj) upper triangle
      for (std::size_t cidx = 0; cidx < work.size(); ++cidx) {
        const auto &w = work[cidx];
        const std::size_t k = w.g.size();
        if (k == 0) continue;

        for (std::size_t ii = 0; ii < k; ++ii) {
          const arma::uword mi = w.g[ii];

          const double h1i_p1 = Hap1_mat(P1, mi);
          const double h2i_p1 = Hap2_mat(P1, mi);
          const double h1i_p2 = Hap1_mat(P2, mi);
          const double h2i_p2 = Hap2_mat(P2, mi);

          for (std::size_t jj = ii; jj < k; ++jj) {
            const arma::uword mj = w.g[jj];
            const double mult = (ii == jj) ? 1.0 : 2.0;

            const double rcf = RC(mi, mj);

            const double h1j_p1 = Hap1_mat(P1, mj);
            const double h2j_p1 = Hap2_mat(P1, mj);
            const double h1j_p2 = Hap1_mat(P2, mj);
            const double h2j_p2 = Hap2_mat(P2, mj);

            auto offdiag_D = [](double a_i, double a_j, double b_i, double b_j) -> double {
              const double mean_i = 0.5 * (a_i + b_i);
              const double mean_j = 0.5 * (a_j + b_j);
              return 0.5 * (a_i * a_j + b_i * b_j) - (mean_i * mean_j);
            };

            const double d1 = offdiag_D(h1i_p1, h1j_p1, h2i_p1, h2j_p1);
            const double d2 = offdiag_D(h1i_p2, h1j_p2, h2i_p2, h2j_p2);

            const double DGen  = rcf * (d1 + d2);
            const double DGen2 = DGen * DGen;

            for (arma::uword ti = 0; ti < numTrait; ++ti) {
              const double UiA = U_mat(mi, ti);
              const double UiD = D_mat(mi, ti);
              for (arma::uword tj = ti; tj < numTrait; ++tj) {
                const arma::uword kcol = tri_u_idx_incl(ti, tj);

                // match your original orientation:
                // U(mj,tj) * DGen * U(mi,ti)
                results1A(x, kcol) += mult * (U_mat(mj, tj) * DGen  * UiA);
                results1D(x, kcol) += mult * (D_mat(mj, tj) * DGen2 * UiD);
              }
            }
          }
        }
      }

      // diagonal outputs
      for (arma::uword ti = 0; ti < numTrait; ++ti) {
        const arma::uword kdiag = tri_u_idx_incl(ti, ti);
        const double varA = results1A(x, kdiag);
        const double varD = results1D(x, kdiag);

        const double eG  = results2(x, OFF_EG  + ti);
        const double eTG = results2(x, OFF_ETG + ti);

        results2(x, OFF_VAR_A + ti) = varA;
        results2(x, OFF_VAR_D + ti) = varD;
        results2(x, OFF_SPV   + ti) = eG  + intensity * std::sqrt(varA);
        results2(x, OFF_TSPV  + ti) = eTG + (std::sqrt(varA) + std::sqrt(varD));
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

      if (calcindex) {
        arma::mat VA = GA;
        arma::mat VD = GD;
        arma::mat L;
        bool spd = arma::chol(L, VA);
        if (!spd) {
          any_psd.store(true, std::memory_order_relaxed);
          arma::vec eval;
          arma::mat evec;
          arma::eig_sym(eval, evec, VA);

          eval.transform([](double x){ return (x < 0.0) ? 0.0 : x; });

          VA = evec * arma::diagmat(eval) * evec.t();
        }

        spd = arma::chol(L, VD);
        if (!spd) {
          any_psd.store(true, std::memory_order_relaxed);
          arma::vec eval;
          arma::mat evec;
          arma::eig_sym(eval, evec, VD);

          eval.transform([](double x){ return (x < 0.0) ? 0.0 : x; });

          VD = evec * arma::diagmat(eval) * evec.t();
        }

        double var_index_A = arma::as_scalar(weights_vec.t() * VA * weights_vec);
        var_index_A = std::max(0.0, var_index_A); // numeric safety
        double index_A = arma::as_scalar(results2.row(x).cols(0, numTrait - 1) * weights_vec);
        double spvi_A  = index_A + intensity * std::sqrt(var_index_A);

        results2(x, OFF_IDX_A)    = index_A;
        results2(x, OFF_VARIDX_A) = var_index_A;
        results2(x, OFF_SPVI_A)   = spvi_A;

        arma::mat VAD = VA + VD;

        double var_index_AD = arma::as_scalar(weights_vec.t() * VAD * weights_vec);
        var_index_AD = std::max(0.0, var_index_AD);
        double index_AD = arma::as_scalar(results2.row(x).cols(numTrait, 2 * numTrait - 1) * weights_vec);
        double spvi_AD  = index_AD + intensity * std::sqrt(var_index_AD);

        results2(x, OFF_IDX_AD)    = index_AD;
        results2(x, OFF_VARIDX_AD) = var_index_AD-var_index_A;
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
      Rcpp::Named("covD")         = out_D,
      Rcpp::Named("check_psd")  = any_psd.load(std::memory_order_relaxed)
    );
  }

  return Rcpp::wrap(results2);
}
