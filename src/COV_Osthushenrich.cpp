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

// Inline function to calculate recombination fraction
inline double qjk(double x, double y, double t) {
  double diff = std::abs(x - y);
  double rcf  = 0.5 * (1.0 - std::exp(-2.0 * diff));
  double oneMinusRcf = 1.0 - rcf;
  return 0.5 + ((1.0 - 2.0 * rcf) / 2.0) * std::pow(oneMinusRcf, t);
}

// [[Rcpp::export]]
SEXP cpp_calculate_covariance_osthushenrich(const NumericMatrix& Crosses,
                  const List& genMap,
                  const NumericMatrix& M,
                  const NumericMatrix& U,
                  double t,
                  double intensity,
                  bool covariance = false,
                  int nThreads = 4) {
#ifdef _OPENMP
  omp_set_num_threads(nThreads);  // set OpenMP threads
#endif

  const arma::uword numCrosses   = Crosses.nrow();
  const arma::uword numMarkers   = M.ncol();
  const arma::uword numTrait     = U.ncol();
  const arma::uword numTraitComb = numTrait * (numTrait + 1) / 2;

  // Convert inputs to Armadillo
  arma::mat M_mat = as<arma::mat>(M);   // (n_individuals × numMarkers)
  arma::mat U_mat = as<arma::mat>(U);   // (numMarkers × numTrait)

  const arma::uword OFF_EG  = 0;
  const arma::uword OFF_VAR = numTrait;
  const arma::uword OFF_SPV = 2 * numTrait;

  // Precompute recombination fractions for each chromosome
  arma::mat QJK(numMarkers, numMarkers, fill::zeros);
  arma::uword startIdx = 0;
  std::vector<std::pair<arma::uword, arma::uword>> chromosomeRanges;

  for (R_xlen_t c = 0; c < genMap.size(); ++c) {
    NumericMatrix gm = as<NumericMatrix>(genMap[c]);  // col 0 = position
    const arma::uword n = gm.nrow();

    for (arma::uword j = 0; j < n; ++j) {
      for (arma::uword k = j; k < n; ++k) {  // upper triangle
        double value = 0.5 * qjk(gm(j, 0), gm(k, 0), t) - 0.25;
        QJK(startIdx + j, startIdx + k) = value;
        QJK(startIdx + k, startIdx + j) = value;  // symmetry
      }
    }
    chromosomeRanges.emplace_back(startIdx, startIdx + n - 1);
    startIdx += n;
  }

  // Precompute M .* mu per trait: P_mu(n_individuals × numMarkers × numTrait)
  const arma::uword nInd = M_mat.n_rows;
  arma::cube P_mu(nInd, numMarkers, numTrait);
  for (arma::uword trait = 0; trait < numTrait; ++trait) {
    arma::rowvec mu_t = U_mat.col(trait).t();        // 1 × numMarkers
    P_mu.slice(trait) = M_mat.each_row() % mu_t;     // n_ind × numMarkers
  }

  // Copy parent indices out of Crosses (1-based in R → 0-based here)
  arma::Col<int> P1_idx(numCrosses), P2_idx(numCrosses);
  for (arma::uword x = 0; x < numCrosses; ++x) {
    P1_idx[x] = static_cast<int>(Crosses(x, 0)) - 1;
    P2_idx[x] = static_cast<int>(Crosses(x, 1)) - 1;
  }

  // Results
  arma::mat results1(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results2(numCrosses, numTrait*3, arma::fill::zeros);

  // Helper to map (ti, tj) with 0 ≤ ti ≤ tj < numTrait to column index
  auto tri_u_idx_incl = [numTrait](arma::uword ti, arma::uword tj) -> arma::uword {
    return ti * numTrait - (ti * (ti - 1)) / 2 + (tj - ti);
  };

#pragma omp parallel for schedule(dynamic)
  for (R_xlen_t x = 0; x < (R_xlen_t)numCrosses; ++x) {
    const int P1 = P1_idx[x];
    const int P2 = P2_idx[x];
    if (P1 < 0 || P2 < 0 || P1 >= (int)nInd || P2 >= (int)nInd) continue;

    // Markers that differ between the two parents (exact compare; use tolerance if needed)
    arma::uvec differing = find(abs(M_mat.row(P1) - M_mat.row(P2)) > 1e-12);
    if (differing.n_elem == 0) continue;

    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      for (arma::uword tj = ti; tj < numTrait; ++tj) {
        if(!covariance && ti!= tj) continue;

        double SigmaSqP1P2 = 0.0;  // reset per trait pair

        // Accumulate by chromosome
        for (const auto& chr : chromosomeRanges) {
          const arma::uword startC = chr.first;
          const arma::uword endC   = chr.second;

          // differing markers within [startC, endC]
          arma::uvec idx    = find( (differing >= startC) % (differing <= endC) );
          if (idx.n_elem == 0) continue;
          arma::uvec chrDiff = differing.elem(idx);

          const arma::uword mlen = chrDiff.n_elem;
          for (arma::uword a = 0; a < mlen; ++a) {
            const arma::uword mi = chrDiff[a];
            const double di = P_mu(P1, mi, ti) - P_mu(P2, mi, ti);

            for (arma::uword b = a; b < mlen; ++b) {
              const arma::uword mj = chrDiff[b];
              const double dj = P_mu(P1, mj, tj) - P_mu(P2, mj, tj);
              const double contrib = di * QJK(mi, mj) * dj;
              SigmaSqP1P2 += (a == b) ? contrib : 2.0 * contrib;
            }
          }
        }

        const arma::uword k = tri_u_idx_incl(ti, tj);
        results1(x, k) = SigmaSqP1P2;

        if (ti == tj) {
          double G1 = arma::dot(M_mat.row(P1), U_mat.col(ti));
          double G2 = arma::dot(M_mat.row(P2), U_mat.col(ti));
          double eG = 0.5 * (G1 + G2);

          results2(x, OFF_EG  + ti) = eG;
          results2(x, OFF_VAR + ti) = SigmaSqP1P2;
          results2(x, OFF_SPV + ti) = eG + intensity * std::sqrt(SigmaSqP1P2);
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
