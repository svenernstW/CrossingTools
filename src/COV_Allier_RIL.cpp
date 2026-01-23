// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <atomic>
#include "parallel_backend.h"

using namespace Rcpp;
using namespace arma;

inline double cjl(double x, double y, int t) {
  double diff = std::abs(x - y);
  double rcf  = 0.5 * (1.0 - std::exp(-2.0 * diff));
  return ((2.0 * rcf) / (1.0 + 2.0 * rcf)) *
    (1.0 - std::pow(0.5, t) * std::pow(1.0 - 2.0 * rcf, t));
}

// packed upper-triangle indexing
inline size_t tri_index(size_t i, size_t j, size_t n) {
  if (i > j) std::swap(i, j);
  return i * n - (i * (i - 1)) / 2 + (j - i);
}
inline double tri_get(const std::vector<double>& tri, size_t i, size_t j, size_t n) {
  return tri[tri_index(i, j, n)];
}
inline void tri_set(std::vector<double>& tri, size_t i, size_t j, size_t n, double v) {
  tri[tri_index(i, j, n)] = v;
}

// [[Rcpp::export]]
SEXP cpp_calculate_covariance_RIL_allier(
    const NumericMatrix& Crosses,
    const List& genMap,
    const NumericMatrix& M,
    const NumericMatrix& U,
    double t,
    double intensity,
    const NumericVector& weights,
    bool covariance = false,
    bool calcindex = false,
    int nThreads = 4) {

  ct_set_threads(nThreads);
  std::atomic<bool> any_psd(false);

  const arma::uword numCrosses   = Crosses.nrow();
  const arma::uword numMarkers   = M.ncol();
  const arma::uword numTrait     = U.ncol();
  const arma::uword numTraitComb = numTrait * (numTrait + 1) / 2;

  arma::mat M_mat = as<arma::mat>(M);   // (n_individuals × numMarkers)
  arma::mat U_mat = as<arma::mat>(U);   // (numMarkers × numTrait)
  const arma::uword nInd = M_mat.n_rows;
  arma::vec weights_vec = as<arma::vec>(weights);

  const arma::uword OFF_EG  = 0;
  const arma::uword OFF_VAR = numTrait;
  const arma::uword OFF_SPV = 2 * numTrait;

  // --- Precompute chromosome ranges + per-chromosome packed CK caches (global, reused by all crosses) ---
  struct ChrCache {
    arma::uword startC;
    arma::uword endC;
    arma::uword nc;
    std::vector<double> CK1;
    std::vector<double> CK2;
  };

  std::vector<ChrCache> chr;
  chr.reserve(genMap.size());

  {
    arma::uword startIdx = 0;
    const int t_int = static_cast<int>(t);

    for (R_xlen_t ci = 0; ci < genMap.size(); ++ci) {
      NumericMatrix gm = as<NumericMatrix>(genMap[ci]);  // col0 = position
      const arma::uword nc = gm.nrow();

      ChrCache cc;
      cc.startC = startIdx;
      cc.endC   = startIdx + nc - 1;
      cc.nc     = nc;

      std::vector<double> pos(nc);
      for (arma::uword j = 0; j < nc; ++j) pos[j] = gm(j, 0);

      const std::size_t L = static_cast<std::size_t>(nc) * (static_cast<std::size_t>(nc) + 1) / 2;
      cc.CK1.assign(L, 0.0);
      cc.CK2.assign(L, 0.0);

      for (arma::uword i = 0; i < nc; ++i) {
        for (arma::uword j = i; j < nc; ++j) {
          double cj_t = cjl(pos[i], pos[j], t_int);
          double cj_1 = cjl(pos[i], pos[j], 1);

          double ck1 = 1.0 - 2.0 * cj_t - std::pow(0.5 * (1.0 - 2.0 * cj_1), t);
          double ck2 = (1.0 - cj_t) * (1.0 - 2.0 * cj_1);

          tri_set(cc.CK1, i, j, static_cast<std::size_t>(nc), ck1);
          tri_set(cc.CK2, i, j, static_cast<std::size_t>(nc), ck2);
        }
      }

      chr.emplace_back(std::move(cc));
      startIdx += nc;
    }

    if (startIdx != numMarkers) {
      stop("Sum of genMap marker counts (%d) != ncol(M) (%d). Check ordering/alignment.",
           (int)startIdx, (int)numMarkers);
    }
  }

  // Copy parent indices out of Crosses (1-based in R → 0-based here)
  arma::Mat<int> P(numCrosses, 4);
  for (arma::uword x = 0; x < numCrosses; ++x) {
    P(x,0) = static_cast<int>(Crosses(x,0)) - 1;
    P(x,1) = static_cast<int>(Crosses(x,1)) - 1;
    P(x,2) = static_cast<int>(Crosses(x,2)) - 1;
    P(x,3) = static_cast<int>(Crosses(x,3)) - 1;
  }

  // Results
  arma::mat results1(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results2(numCrosses, numTrait * 3 + 3, arma::fill::zeros);

  auto tri_u_idx_incl = [numTrait](arma::uword ti, arma::uword tj) -> arma::uword {
    return ti * numTrait - (ti * (ti - 1)) / 2 + (tj - ti);
  };


  ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
    arma::uword x = static_cast<arma::uword>(xi);

    const int P1 = P(x,0), P2 = P(x,1), P3 = P(x,2), P4 = P(x,3);
    if (P1 < 0 || P2 < 0 || P3 < 0 || P4 < 0 ||
        P1 >= (int)nInd || P2 >= (int)nInd || P3 >= (int)nInd || P4 >= (int)nInd) {
      return;
    }

    // EGBVs once per trait (genome-wide)
    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      double G1 = arma::dot(M_mat.row(P1), U_mat.col(ti));
      double G2 = arma::dot(M_mat.row(P2), U_mat.col(ti));
      double G3 = arma::dot(M_mat.row(P3), U_mat.col(ti));
      double G4 = arma::dot(M_mat.row(P4), U_mat.col(ti));
      results2(x, OFF_EG + ti) = 0.25 * (G1 + G2 + G3 + G4);
    }

    // ------------------------------------------------------------
    // Precompute per-chromosome work + detect if ANY marker differs
    // ------------------------------------------------------------
    bool any_diff = false;

    struct ChrWork {
      std::vector<arma::uword> li;
      std::vector<double> a12, a34, a14, a13, a24, a23;
    };

    std::vector<ChrWork> work(chr.size());

    for (std::size_t cidx = 0; cidx < chr.size(); ++cidx) {
      const auto& cc = chr[cidx];
      const arma::uword startC = cc.startC;
      const arma::uword endC   = cc.endC;

      auto& w = work[cidx];
      w.li.clear();
      w.a12.clear(); w.a34.clear();
      w.a14.clear(); w.a13.clear();
      w.a24.clear(); w.a23.clear();

      w.li.reserve(256);
      w.a12.reserve(256); w.a34.reserve(256);
      w.a14.reserve(256); w.a13.reserve(256);
      w.a24.reserve(256); w.a23.reserve(256);

      for (arma::uword g = startC; g <= endC; ++g) {
        const double m1 = M_mat(P1, g);
        const double m2 = M_mat(P2, g);
        const double m3 = M_mat(P3, g);
        const double m4 = M_mat(P4, g);

        // if all equal, skip
        if (m1 == m2 && m1 == m3 && m1 == m4) continue;

        any_diff = true;

        const arma::uword li = g - startC;
        w.li.push_back(li);

        w.a12.push_back(m1 - m2);
        w.a34.push_back(m3 - m4);
        w.a14.push_back(m1 - m4);
        w.a13.push_back(m1 - m3);
        w.a24.push_back(m2 - m4);
        w.a23.push_back(m2 - m3);
      }
    }

    // ---- NEW: if no differences anywhere, Var=0 and SPV=eG ----
    if (!any_diff) {
      for (arma::uword ti = 0; ti < numTrait; ++ti) {
        const double eG = results2(x, OFF_EG + ti);
        results2(x, OFF_VAR + ti) = 0.0;
        results2(x, OFF_SPV + ti) = eG;
      }
      // results1(x,*) stays 0 => V is 0 later => index variance 0 => spvi=index
      return;
    }

    // ------------------------------------------------------------
    // Trait pairs (reuse work[cidx])
    // ------------------------------------------------------------
    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      for (arma::uword tj = ti; tj < numTrait; ++tj) {
        if (!covariance && ti != tj) continue;

        double Sigma = 0.0;

        for (std::size_t cidx = 0; cidx < chr.size(); ++cidx) {
          const auto& cc = chr[cidx];
          const arma::uword startC = cc.startC;
          const arma::uword nc     = cc.nc;

          const auto& w = work[cidx];
          const std::size_t k = w.li.size();
          if (k == 0) continue;

          double add_chr = 0.0;

          for (std::size_t ii = 0; ii < k; ++ii) {
            const arma::uword li = w.li[ii];
            const arma::uword gi = startC + li;

            const double a12_i = w.a12[ii];
            const double a34_i = w.a34[ii];
            const double a14_i = w.a14[ii];
            const double a13_i = w.a13[ii];
            const double a24_i = w.a24[ii];
            const double a23_i = w.a23[ii];

            for (std::size_t jj = ii; jj < k; ++jj) {
              const arma::uword lj = w.li[jj];
              const arma::uword gj = startC + lj;

              const double D12  = 0.0625 * (a12_i * w.a12[jj]);
              const double D34  = 0.0625 * (a34_i * w.a34[jj]);
              const double phi1 = D12 + D34;

              const double D14  = 0.0625 * (a14_i * w.a14[jj]);
              const double D13  = 0.0625 * (a13_i * w.a13[jj]);
              const double D24  = 0.0625 * (a24_i * w.a24[jj]);
              const double D23  = 0.0625 * (a23_i * w.a23[jj]);
              const double phi2 = D14 + D13 + D24 + D23;

              const double ck1 = tri_get(cc.CK1, (std::size_t)li, (std::size_t)lj, (std::size_t)nc);
              const double ck2 = tri_get(cc.CK2, (std::size_t)li, (std::size_t)lj, (std::size_t)nc);

              const double Dcomb = (ck1 * phi2) + (ck2 * phi1);

              const double contrib = U_mat(gi, ti) * Dcomb * U_mat(gj, tj);
              add_chr += (ii == jj) ? contrib : 2.0 * contrib;
            }
          }

          Sigma += add_chr;
        }

        const arma::uword kcol = tri_u_idx_incl(ti, tj);
        results1(x, kcol) = Sigma;

        if (ti == tj) {
          const double eG = results2(x, OFF_EG + ti);
          results2(x, OFF_VAR + ti) = Sigma;
          results2(x, OFF_SPV + ti) = eG + intensity * std::sqrt(Sigma);
        }
      }
    }
  });



  // If requested, convert each row to an nTrait×nTrait symmetric matrix + desired gains index
  if (covariance) {
    std::vector<arma::mat> covs(numCrosses);

    ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
      arma::uword x = static_cast<arma::uword>(xi);

      arma::mat V(numTrait, numTrait, arma::fill::zeros);
      for (arma::uword ti = 0; ti < numTrait; ++ti) {
        for (arma::uword tj = ti; tj < numTrait; ++tj) {
          const arma::uword k = tri_u_idx_incl(ti, tj);
          const double v = results1(x, k);
          V(ti, tj) = v;
          V(tj, ti) = v;
        }
      }
      covs[x] = V;

      if (calcindex) {
        arma::mat L;
        bool spd = arma::chol(L, V);
        if (!spd) {
          any_psd.store(true, std::memory_order_relaxed);
          arma::vec eval;
          arma::mat evec;
          arma::eig_sym(eval, evec, V);

          eval.transform([](double x){ return (x < 0.0) ? 0.0 : x; });

          V = evec * arma::diagmat(eval) * evec.t();
        }
        double sigma_gains = arma::as_scalar(weights_vec.t() * V * weights_vec);
        double index = arma::as_scalar(results2.row(x).cols(0, numTrait - 1) * weights_vec);
        double spvi  = index + intensity * std::sqrt(sigma_gains);

        results2(x, 3 * numTrait + 0) = index;
        results2(x, 3 * numTrait + 1) = sigma_gains;
        results2(x, 3 * numTrait + 2) = spvi;
      }
    });

    Rcpp::List out(numCrosses);
    for (arma::uword x = 0; x < numCrosses; ++x) out[x] = Rcpp::wrap(covs[x]);

    return Rcpp::List::create(
      Rcpp::Named("cross_values") = results2,
      Rcpp::Named("covariances")  = out,
      Rcpp::Named("check_psd")  = any_psd.load(std::memory_order_relaxed)
    );
  }

  return Rcpp::wrap(results2);
}
