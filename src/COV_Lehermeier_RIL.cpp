// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace arma;

// Inline function to calculate recombination fraction
inline double qjk(double x, double y, double t) {
  double diff = std::abs(x - y);
  double rcf = 0.5 * (1 - std::exp(-2 * diff));
  double ck = ((2 * rcf)/(1 + 2 * rcf)) *((1- std::pow(0.5 ,t))*std::pow(1 - rcf, t));

    return 2 * ck - std::pow(0.5 * (1-2 * rcf),t);
}

// [[Rcpp::export]]
SEXP cpp_calculate_covariance_RIL_lehermeier(const NumericMatrix& Crosses,
                                             const List& genMap,
                                             const NumericMatrix& M,
                                             const NumericMatrix& U,
                                             double t,
                                             double intensity,
                                             const NumericVector& gains,
                                             bool covariance = false,
                                             bool calcgains = false,
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
  const arma::uword nInd = M_mat.n_rows;
  arma::vec gains_vec = as<arma::vec>(gains);  // length == numTrait

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
  if (startIdx != numMarkers)
    stop("Sum of genMap marker counts (%d) != numMarkers (%d).", (int)startIdx, (int)numMarkers);


  // Copy parent indices out of Crosses (1-based in R → 0-based here)
  arma::Col<int> P1_idx(numCrosses), P2_idx(numCrosses);
  for (arma::uword x = 0; x < numCrosses; ++x) {
    P1_idx[x] = static_cast<int>(Crosses(x, 0)) - 1;
    P2_idx[x] = static_cast<int>(Crosses(x, 1)) - 1;
  }

  // Results
  arma::mat results1(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results2(numCrosses, numTrait*3+3, arma::fill::zeros);

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

        double SigmaSqP1P2 = 0.0;

        // Iterate through each chromosome range
        for (const auto& chrRange : chromosomeRanges) {
          int startC = chrRange.first;
          int endC   = chrRange.second;

          // differing markers within [startC, endC]
          arma::uvec idx    = find( (differing >= startC) % (differing <= endC) );
          if (idx.n_elem == 0) continue;
          arma::uvec chrDiff = differing.elem(idx);

          // If no differing markers on this chromosome, continue
          if (chrDiff.n_elem == 0) {
            continue;
          }

          // Now only loop over these chromosome-specific markers
          for (uword i = 0; i < chrDiff.n_elem; ++i) {
            int marker_i = chrDiff[i];
            for (uword j = i; j < chrDiff.n_elem; ++j) {
              int marker_j = chrDiff[j];
              double Dprime = (0.0625 * ((M_mat(P1, marker_i)-M_mat(P2, marker_i))*(M_mat(P1, marker_j)-M_mat(P2, marker_j))));

              double D = (4 * Dprime) * (1 -  QJK(marker_i, marker_j)); //Note that different from DH function the times 2 is already in the QJK function

              double contrib = U_mat(marker_i, ti) * D * U_mat(marker_j, tj);

              SigmaSqP1P2 += (i == j) ? contrib : 2 * contrib;
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
    std::vector<arma::mat> covs(numCrosses);

#pragma omp parallel for schedule(dynamic)
    for (arma::uword x = 0; x < numCrosses; ++x) {
      arma::mat G(numTrait, numTrait, arma::fill::zeros);
      for (arma::uword ti = 0; ti < numTrait; ++ti) {
        for (arma::uword tj = ti; tj < numTrait; ++tj) {
          const arma::uword k = tri_u_idx_incl(ti, tj); // 0..Ktri-1
          const double v = results1(x, k);
          G(ti, tj) = v;
          G(tj, ti) = v;
        }
      }
      covs[x] = std::move(G);

      if (calcgains) {
        arma::mat V = covs[x];
        arma::mat L;
        bool spd = arma::chol(L, V);
        if (!spd) {
          arma::vec eval; arma::mat evec;
          arma::eig_sym(eval, evec, V);                 // O(n^3) but n is tiny
          double tol = std::max(1e-12, 1e-8 * eval.max());
          for (auto& l : eval) if (l < tol) l = tol;    // clamp
          V = evec * arma::diagmat(eval) * evec.t();    // PSD → SPD (since tol>0)
          arma::chol(L, V);                             // must succeed now
        }


        arma::vec w;
        arma::solve(w, V, gains_vec, arma::solve_opts::likely_sympd);

        double sigma_gains = arma::as_scalar(w.t() * V * w);

        double index =
          arma::as_scalar( results2.row(x).cols(0, numTrait - 1) * w );

        double spvi = index + intensity * std::sqrt(sigma_gains);

        results2(x, 3 * numTrait + 0) = index;
        results2(x, 3 * numTrait + 1) = sigma_gains;
        results2(x, 3 * numTrait + 2) = spvi;
      }
    }

    // wrap covariances into an R list
    Rcpp::List out(numCrosses);
    for (arma::uword x = 0; x < numCrosses; ++x) out[x] = Rcpp::wrap(covs[x]);

    return Rcpp::List::create(
      Rcpp::Named("cross_values") = results2,
      Rcpp::Named("covariances")  = out
    );
  }

  return Rcpp::wrap(results2);    // returns upper triangle of the covariance matrix in vectorized form
}
