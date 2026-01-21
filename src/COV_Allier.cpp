// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include "parallel_backend.h"

using namespace Rcpp;
using namespace arma;

// Recombination helper
inline double cjl(double x, double y, int t) {
  double diff = std::abs(x - y);
  double rcf  = 0.5 * (1.0 - std::exp(-2.0 * diff));
  return ((2.0 * rcf) / (1.0 + 2.0 * rcf)) *
    (1.0 - std::pow(0.5, t) * std::pow(1.0 - 2.0 * rcf, t));
}

// Packed upper-triangle indexing for size n
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
SEXP cpp_calculate_covariance_allier(
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

  const arma::uword numCrosses   = Crosses.nrow();
  const arma::uword numMarkers   = M.ncol();
  const arma::uword numTrait     = U.ncol();
  const arma::uword numTraitComb = numTrait * (numTrait + 1) / 2;

  arma::mat M_mat = as<arma::mat>(M);   // (n_individuals × numMarkers)
  arma::mat U_mat = as<arma::mat>(U);   // (numMarkers × numTrait)
  const arma::uword nInd = M_mat.n_rows;
  arma::vec weights_vec = as<arma::vec>(weights); // length == numTrait

  const arma::uword OFF_EG  = 0;
  const arma::uword OFF_VAR = numTrait;
  const arma::uword OFF_SPV = 2 * numTrait;

  // --- Precompute chromosome ranges + per-chromosome packed recombination caches ---
  struct ChrCache {
    arma::uword startC;
    arma::uword endC;
    arma::uword nc;
    std::vector<double> CK;
    std::vector<double> CK1;
    std::vector<double> C1;
  };

  std::vector<ChrCache> chr;
  chr.reserve(genMap.size());

  {
    arma::uword startIdx = 0;
    const int t_int = static_cast<int>(t);
    const int tm1   = t_int - 1;

    for (R_xlen_t ci = 0; ci < genMap.size(); ++ci) {
      NumericMatrix gm = as<NumericMatrix>(genMap[ci]);  // col0 = pos
      const arma::uword nc = gm.nrow();

      ChrCache cc;
      cc.startC = startIdx;
      cc.endC   = startIdx + nc - 1;
      cc.nc     = nc;

      const std::size_t L = static_cast<std::size_t>(nc) * (static_cast<std::size_t>(nc) + 1) / 2;
      cc.CK.assign(L,  0.0);
      cc.CK1.assign(L, 0.0);
      cc.C1.assign(L,  0.0);

      // positions
      std::vector<double> pos(nc);
      for (arma::uword j = 0; j < nc; ++j) pos[j] = gm(j, 0);

      for (arma::uword i = 0; i < nc; ++i) {
        for (arma::uword j = i; j < nc; ++j) {
          double cj_t   = cjl(pos[i], pos[j], t_int);
          double cj_tm1 = cjl(pos[i], pos[j], tm1);
          double cj_1   = cjl(pos[i], pos[j], 1);

          double ck  = 1.0 - 2.0 * cj_t;   // CK
          double ck1 = cj_tm1;             // CK1
          double c1  = 1.0 - 2.0 * cj_1;   // C1

          tri_set(cc.CK,  i, j, static_cast<std::size_t>(nc), ck);
          tri_set(cc.CK1, i, j, static_cast<std::size_t>(nc), ck1);
          tri_set(cc.C1,  i, j, static_cast<std::size_t>(nc), c1);
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

    // ---- (A) Precompute eG once per trait (same as your first code) ----
    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      double G1 = arma::dot(M_mat.row(P1), U_mat.col(ti));
      double G2 = arma::dot(M_mat.row(P2), U_mat.col(ti));
      double G3 = arma::dot(M_mat.row(P3), U_mat.col(ti));
      double G4 = arma::dot(M_mat.row(P4), U_mat.col(ti));
      results2(x, OFF_EG + ti) = 0.25 * (G1 + G2 + G3 + G4);
    }

    // ---- (B) Precompute diff markers ONCE per chromosome for this cross ----
    std::vector<std::vector<arma::uword>> diff_by_chr;
    diff_by_chr.resize(chr.size());

    for (std::size_t cidx = 0; cidx < chr.size(); ++cidx) {
      const auto& cc = chr[cidx];
      diff_by_chr[cidx].clear();
      diff_by_chr[cidx].reserve(256);

      for (arma::uword g = cc.startC; g <= cc.endC; ++g) {
        double a = M_mat(P1, g), b = M_mat(P2, g), c = M_mat(P3, g), d = M_mat(P4, g);
        // keep marker if not all equal
        if ((a!=b) || (a!=c) || (a!=d) || (b!=c) || (b!=d) || (c!=d)) {
          diff_by_chr[cidx].push_back(g);
        }
      }
    }

    // ---- (C) Now trait pairs reuse the precomputed diff lists ----
    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      for (arma::uword tj = ti; tj < numTrait; ++tj) {
        if (!covariance && ti != tj) continue;

        double Sigma = 0.0;

        for (std::size_t cidx = 0; cidx < chr.size(); ++cidx) {
          const auto& cc   = chr[cidx];
          const auto& diff = diff_by_chr[cidx];
          if (diff.empty()) continue;

          const arma::uword startC = cc.startC;
          const arma::uword nc     = cc.nc;

          double add_chr = 0.0;

          for (std::size_t ii = 0; ii < diff.size(); ++ii) {
            arma::uword gi = diff[ii];
            arma::uword li = gi - startC;

            for (std::size_t jj = ii; jj < diff.size(); ++jj) {
              arma::uword gj = diff[jj];
              arma::uword lj = gj - startC;

              double D12 = 0.0625 * ((M_mat(P1,gi)-M_mat(P2,gi))*(M_mat(P1,gj)-M_mat(P2,gj)));
              double D34 = 0.0625 * ((M_mat(P3,gi)-M_mat(P4,gi))*(M_mat(P3,gj)-M_mat(P4,gj)));
              double phi1 = D12 + D34;

              double D14 = 0.0625 * ((M_mat(P1,gi)-M_mat(P4,gi))*(M_mat(P1,gj)-M_mat(P4,gj)));
              double D13 = 0.0625 * ((M_mat(P1,gi)-M_mat(P3,gi))*(M_mat(P1,gj)-M_mat(P3,gj)));
              double D24 = 0.0625 * ((M_mat(P2,gi)-M_mat(P4,gi))*(M_mat(P2,gj)-M_mat(P4,gj)));
              double D23 = 0.0625 * ((M_mat(P2,gi)-M_mat(P3,gi))*(M_mat(P2,gj)-M_mat(P3,gj)));
              double phi2 = D14 + D13 + D24 + D23;

              double ck  = tri_get(cc.CK,  (std::size_t)li, (std::size_t)lj, (std::size_t)nc);
              double ck1 = tri_get(cc.CK1, (std::size_t)li, (std::size_t)lj, (std::size_t)nc);
              double c1  = tri_get(cc.C1,  (std::size_t)li, (std::size_t)lj, (std::size_t)nc);

              double Dcomb   = (ck * phi2) + ((ck + ck1) * c1 * phi1);
              double contrib = U_mat(gi, ti) * Dcomb * U_mat(gj, tj);

              add_chr += (ii == jj) ? contrib : 2.0 * contrib;
            }
          }

          Sigma += add_chr;
        }

        const arma::uword k = tri_u_idx_incl(ti, tj);
        results1(x, k) = Sigma;

        if (ti == tj) {
          double eG = results2(x, OFF_EG + ti);
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
      Rcpp::Named("covariances")  = out
    );
  }

  return Rcpp::wrap(results2);
}
