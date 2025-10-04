// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// ---------------------------------------------
// Recombination helper (same math as before)
// ---------------------------------------------
inline double cjl(double x, double y, int t) {
  double diff = std::abs(x - y);
  double rcf  = 0.5 * (1.0 - std::exp(-2.0 * diff));
  return ((2.0 * rcf) / (1.0 + 2.0 * rcf)) *
    (1.0 - std::pow(0.5, t) * std::pow(1.0 - 2.0 * rcf, t));
}

// ---------------------------------------------
// Packed upper-triangle indexing for size n
// Stores entries for i<=j in a flat vector:
// idx(i,j) = number of elements before row i +
//            offset (j - i) in that row
// Elements before row i = i*n - i*(i-1)/2
// ---------------------------------------------
inline size_t tri_index(size_t i, size_t j, size_t n) {
  if (i > j) std::swap(i, j);
  // i <= j now
  return i * n - (i * (i - 1)) / 2 + (j - i);
}
inline double tri_get(const std::vector<double>& tri, size_t i, size_t j, size_t n) {
  return tri[tri_index(i, j, n)];
}
inline void tri_set(std::vector<double>& tri, size_t i, size_t j, size_t n, double v) {
  tri[tri_index(i, j, n)] = v;
}

// [[Rcpp::export]]
SEXP cpp_calculate_covariance_allier(const NumericMatrix& Crosses,
             const List& genMap,
             const NumericMatrix& M,
             const NumericMatrix& U,
             double t,
             double intensity,
             bool covariance = false,
             int nThreads = 4) {
#ifdef _OPENMP
  omp_set_dynamic(1);
  omp_set_num_threads(nThreads);
#endif

  const arma::uword numCrosses   = Crosses.nrow();
  const arma::uword numTrait     = U.ncol();
  const arma::uword numTraitComb = numTrait * (numTrait + 1) / 2;

  arma::mat M_mat = as<arma::mat>(M);   // (n_individuals × numMarkers)
  arma::mat U_mat = as<arma::mat>(U);   // (numMarkers × numTrait)

  const arma::uword OFF_EG  = 0;
  const arma::uword OFF_VAR = numTrait;
  const arma::uword OFF_SPV = 2 * numTrait;

  // --- Precompute chromosome ranges + positions ---
  std::vector<std::pair<int,int>> chrRanges;
  std::vector<std::vector<double>> chrPos;
  {
    int startIdx = 0;
    chrRanges.reserve(genMap.size());
    chrPos.reserve(genMap.size());
    for (int i = 0; i < genMap.size(); ++i) {
      NumericMatrix gm = as<NumericMatrix>(genMap[i]);
      int n = gm.nrow();
      chrRanges.emplace_back(startIdx, startIdx + n - 1);
      std::vector<double> pos(n);
      for (int j = 0; j < n; ++j) pos[j] = gm(j, 0);
      chrPos.emplace_back(std::move(pos));
      startIdx += n;
    }
  }

  // Results
  arma::mat results1(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results2(numCrosses, numTrait*3, arma::fill::zeros);

  // Helper to map (ti, tj) with 0 ≤ ti ≤ tj < numTrait to column index
  auto tri_u_idx_incl = [numTrait](arma::uword ti, arma::uword tj) -> arma::uword {
    return ti * numTrait - (ti * (ti - 1)) / 2 + (tj - ti);
  };

  // --- Stream over chromosomes ---
  for (std::size_t chrIdx = 0; chrIdx < chrRanges.size(); ++chrIdx) {
    const int startC = chrRanges[chrIdx].first;
    const int endC   = chrRanges[chrIdx].second;
    const int nc     = endC - startC + 1;

    // Packed arrays length = nc*(nc+1)/2
    const std::size_t L = static_cast<std::size_t>(nc) * (static_cast<std::size_t>(nc) + 1) / 2;
    std::vector<double> CK_pack(L), CK1_pack(L), C1_pack(L);

    const int tm1 = static_cast<int>(t) - 1;
    for (int i = 0; i < nc; ++i) {
      for (int j = i; j < nc; ++j) {
        double pi = chrPos[chrIdx][i];
        double pj = chrPos[chrIdx][j];

        double cj_t   = cjl(pi, pj, static_cast<int>(t));
        double cj_tm1 = cjl(pi, pj, tm1);
        double cj_1   = cjl(pi, pj, 1);

        double ck  = 1.0 - 2.0 * cj_t;   // CK
        double ck1 = cj_tm1;             // CK1
        double c1  = 1.0 - 2.0 * cj_1;   // C1

        tri_set(CK_pack,  i, j, nc, ck);
        tri_set(CK1_pack, i, j, nc, ck1);
        tri_set(C1_pack,  i, j, nc, c1);
      }
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (arma::sword x = 0; x < static_cast<arma::sword>(numCrosses); ++x) {
      int P1 = Crosses(x,0)-1, P2 = Crosses(x,1)-1, P3 = Crosses(x,2)-1, P4 = Crosses(x,3)-1;

      // bounds check against M_mat.n_rows (arma::uword)
      arma::uword nrows = M_mat.n_rows;
      if (P1 < 0 || P2 < 0 || P3 < 0 || P4 < 0 ||
          static_cast<arma::uword>(P1) >= nrows ||
          static_cast<arma::uword>(P2) >= nrows ||
          static_cast<arma::uword>(P3) >= nrows ||
          static_cast<arma::uword>(P4) >= nrows) {
        continue;
      }

      const rowvec m1 = M_mat.row(P1);
      const rowvec m2 = M_mat.row(P2);
      const rowvec m3 = M_mat.row(P3);
      const rowvec m4 = M_mat.row(P4);

      std::vector<int> chrDiff;
      chrDiff.reserve(256);
      for (int g = startC; g <= endC; ++g) {
        double a = m1[g], b = m2[g], c = m3[g], d = m4[g];
        if ((a!=b) || (a!=c) || (a!=d) || (b!=c) || (b!=d) || (c!=d)) chrDiff.push_back(g);
      }
      if (chrDiff.empty()) continue;

      for (arma::uword ti = 0; ti < numTrait; ++ti) {
        for (arma::uword tj = ti; tj < numTrait; ++tj) {
          if(!covariance && ti!= tj) continue;

          double add = 0.0;
          for (std::size_t ii = 0; ii < chrDiff.size(); ++ii) {
            const int gi = chrDiff[ii];
            const int li = gi - startC; // local 0..nc-1

            for (std::size_t jj = ii; jj < chrDiff.size(); ++jj) {
              const int gj = chrDiff[jj];
              const int lj = gj - startC;

              double D12 = 0.0625 * ((M_mat(P1,gi)-M_mat(P2,gi))*(M_mat(P1,gj)-M_mat(P2,gj)));
              double D34 = 0.0625 * ((M_mat(P3,gi)-M_mat(P4,gi))*(M_mat(P3,gj)-M_mat(P4,gj)));
              double phi1 = D12 + D34;

              double D14 = 0.0625 * ((M_mat(P1,gi)-M_mat(P4,gi))*(M_mat(P1,gj)-M_mat(P4,gj)));
              double D13 = 0.0625 * ((M_mat(P1,gi)-M_mat(P3,gi))*(M_mat(P1,gj)-M_mat(P3,gj)));
              double D24 = 0.0625 * ((M_mat(P2,gi)-M_mat(P4,gi))*(M_mat(P2,gj)-M_mat(P4,gj)));
              double D23 = 0.0625 * ((M_mat(P2,gi)-M_mat(P3,gi))*(M_mat(P2,gj)-M_mat(P3,gj)));
              double phi2 = D14 + D13 + D24 + D23;

              double ck  = tri_get(CK_pack,  static_cast<std::size_t>(li), static_cast<std::size_t>(lj), static_cast<std::size_t>(nc));
              double ck1 = tri_get(CK1_pack, static_cast<std::size_t>(li), static_cast<std::size_t>(lj), static_cast<std::size_t>(nc));
              double c1  = tri_get(C1_pack,  static_cast<std::size_t>(li), static_cast<std::size_t>(lj), static_cast<std::size_t>(nc));

              double Dcomb   = (ck * phi2) + ((ck + ck1) * c1 * phi1);

              // Use marker indices gi/gj for U
              double contrib = U_mat(gi, ti) * Dcomb * U_mat(gj, tj);
              add += (ii == jj) ? contrib : 2.0 * contrib;
            }
          }

          const arma::uword k = tri_u_idx_incl(ti, tj);
          results1(x, k) = add;

          if (ti == tj) {
            double G1 = arma::dot(M_mat.row(P1), U_mat.col(ti));
            double G2 = arma::dot(M_mat.row(P2), U_mat.col(ti));
            double G3 = arma::dot(M_mat.row(P3), U_mat.col(ti));
            double G4 = arma::dot(M_mat.row(P4), U_mat.col(ti));

            double eG = 0.25 * (G1 + G2 + G3 + G4);

            results2(x, OFF_EG  + ti) = eG;
            results2(x, OFF_VAR + ti) = add;
            results2(x, OFF_SPV + ti) = eG + intensity * std::sqrt(add);
          }
        }
      }
    }
  }
    // If requested, convert each row to an nTrait×nTrait symmetric matrix
    if (covariance) {
      List out(numCrosses);
      for (uword x = 0; x < numCrosses; ++x) {
        NumericMatrix G(numTrait, numTrait);
        for (uword ti = 0; ti < numTrait; ++ti) {
          for (uword tj = ti; tj < numTrait; ++tj) {
            const uword k = tri_u_idx_incl(ti, tj);
            double v = results1(x, k);
            G(ti, tj) = v;
            G(tj, ti) = v;  // mirror
          }
        }
        out[x] = G;
      }
      return List::create(
        Named("cross_values") = results2,
        Named("covariances") = out
      );
    }

    return Rcpp::wrap(results2);    // returns upper triangle of the covariance matrix in vectorized form
  }
