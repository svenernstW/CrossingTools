// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <atomic>
#include <algorithm>
#include <limits>
#include <utility>
#include "parallel_backend.h"  // OpenMP when available, else RcppParallel

using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
SEXP cpp_calculate_covariance_allier_ril(const NumericMatrix& Crosses,
                   const List& genMap,
                   const NumericMatrix& M,
                   const NumericMatrix& U,
                   int t,
                   int s,
                   double intensity,
                   const NumericVector& weights,
                   const NumericVector& p,
                   bool covariance = false,
                   bool calcindex = false,
                   int nThreads = 4) {

  ct_set_threads(nThreads);
  std::atomic<bool> any_psd(false);

  const arma::uword numCrosses   = Crosses.nrow();
  const arma::uword numMarkers   = M.ncol();
  const arma::uword numTrait     = U.ncol();
  const arma::uword numTraitComb = numTrait * (numTrait + 1) / 2;

  if (Crosses.ncol() < 4 || U.nrow() != numMarkers ||
      p.size() != static_cast<R_xlen_t>(numMarkers) ||
      weights.size() != static_cast<R_xlen_t>(numTrait)) {
    Rcpp::stop("Invalid Crosses, M, U, p, or weights dimensions");
  }
  if (t < 0 || s < 0 || t > std::numeric_limits<int>::max() - s - 2) {
    Rcpp::stop("t and s must be non-negative and their sum must be representable");
  }

  // Convert inputs to Armadillo
  arma::mat M_mat = as<arma::mat>(M);  // genotypes coded 0..2
  arma::mat U_mat = as<arma::mat>(U);  // additive marker effects
  const arma::uword nInd = M_mat.n_rows;
  arma::vec weights_vec = as<arma::vec>(weights);
  arma::rowvec p_vec = as<arma::rowvec>(p);

  // Centred genotype matrix and parental breeding values (unchanged)
  arma::mat W_mat = M_mat;
  W_mat.each_row() -= 2.0 * p_vec;
  arma::mat GEBV = W_mat * U_mat;

  const arma::uword OFF_EG  = 0;
  const arma::uword OFF_VAR = numTrait;
  const arma::uword OFF_SPV = 2 * numTrait;

  // Ordered finite-generation RIL model: t random-mating rounds from the
  // four-way pseudo-F1, followed by s rounds of selfing.
  // At t=0, s=g, this reproduces Allier's finite RIL(g) coefficients.

  // Precompute the TWO four-way Haldane kernels as exponential weights.
  // Rather than multiply a random-mating formula by a selfing formula,
  // advance the two-locus transmissions IN THE SPECIFIED ORDER.
  //
  // Let f be the association between loci on the same homolog, and
  //     h be the association between loci on opposite homologs.
  // One meiosis gives ((1-r)*f + r*h).
  // With z=1-2*r, (1-r)=(1+z)/2 and r=(1-z)/2:
  //   gamete(z) = 0.5*(f+h) + 0.5*z*(f-h).
  //
  // Random mating: f_new = gamete, h_new = 0.
  // Selfing:       f_new = gamete, h_new = 0.5*(f+h).
  //
  // Across blocks (B): fB=1, hB=-1 initially (two fixed F1 homologs).
  // Within blocks (W): fW=z, hW=0 initially (two independently
  //   recombined F1 parental gametes).
  //
  // At the final generation:
  //   cB(z) = (fB(z) + hB(z)) / 8
  //   cW(z) = (fW(z) + hW(z)) / 4
  // These describe diploid RIL-generation progeny, not DH lines.
  // The extra factor 0.25 in coeffB/coeffW is the OUTER factor
  // in the four-way genotype covariance equation from the appendix.

  const int K = t + s + 1;
  std::vector<double> coeffB(K, 0.0), coeffW(K, 0.0);
  double diagonalB = 0.0, diagonalW = 0.0;

  std::vector<double> fB(K+1, 0.0), hB(K+1, 0.0);
  std::vector<double> fW(K+1, 0.0), hW(K+1, 0.0);
  fB[0] = 1.0;
  hB[0] = -1.0;
  fW[1] = 1.0;

  for (int generation = 0; generation < t+s; ++generation) {
    std::vector<double> newFB(K+1, 0.0), newHB(K+1, 0.0);
    std::vector<double> newFW(K+1, 0.0), newHW(K+1, 0.0);
    for (int k = 0; k <= K; ++k) {
      newFB[k] = 0.5*(fB[k]+hB[k]);
      newFW[k] = 0.5*(fW[k]+hW[k]);
      if (k > 0) {
        newFB[k] += 0.5*(fB[k-1]-hB[k-1]);
        newFW[k] += 0.5*(fW[k-1]-hW[k-1]);
      }
      if (generation >= t) {  // selfing AFTER random mating
        newHB[k] = 0.5*(fB[k]+hB[k]);
        newHW[k] = 0.5*(fW[k]+hW[k]);
      }
    }
    fB.swap(newFB); hB.swap(newHB);
    fW.swap(newFW); hW.swap(newHW);
  }

  // Include the outer 0.25 in the four-way covariance equation.
  for (int k = 1; k <= K; ++k) {
    coeffB[k-1] = 0.25 * 0.125 * (fB[k]+hB[k]);
    coeffW[k-1] = 0.25 * 0.250 * (fW[k]+hW[k]);
  }
  for (int k = 0; k < K; ++k) {
    diagonalB += coeffB[k];
    diagonalW += coeffW[k];
  }

  // Precompute marker positions and chromosome ranges (no QJK or P_mu).
  arma::vec markerPos(numMarkers);
  std::vector<std::pair<arma::uword, arma::uword>> chromosomeRanges;
  arma::uword startIdx = 0;

  for (R_xlen_t c = 0; c < genMap.size(); ++c) {
    NumericMatrix gm = as<NumericMatrix>(genMap[c]);
    const arma::uword n = gm.nrow();
    if (n == 0) continue;
    if (gm.ncol() < 1 || startIdx + n > numMarkers) {
      Rcpp::stop("Invalid genetic map dimensions");
    }
    for (arma::uword j = 0; j < n; ++j) {
      if (!std::isfinite(gm(j, 0)) ||
          (j > 0 && gm(j, 0) < gm(j - 1, 0))) {
        Rcpp::stop("Marker positions must be finite and sorted within chromosomes");
      }
      markerPos[startIdx + j] = gm(j, 0); // positions in Morgans
    }
    chromosomeRanges.emplace_back(startIdx, startIdx + n - 1);
    startIdx += n;
  }
  if (startIdx != numMarkers) {
    Rcpp::stop("genMap and M contain different marker counts");
  }

  // Parent indices (1-based in R, 0-based in C++)
  arma::Mat<int> P(numCrosses, 4);
  for (arma::uword x = 0; x < numCrosses; ++x) {
    for (arma::uword a = 0; a < 4; ++a) {
      P(x,a) = static_cast<int>(Crosses(x,a)) - 1;
    }
  }

  // Results: unchanged format
  arma::mat results1(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results2(numCrosses, numTrait * 3 + 3, arma::fill::zeros);

  auto tri_u_idx_incl = [numTrait](arma::uword ti, arma::uword tj) -> arma::uword {
    return ti * numTrait - (ti * (ti - 1)) / 2 + (tj - ti);
  };

  // Parallel over four-way crosses
  ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
    arma::uword x = static_cast<arma::uword>(xi);
    const int P1 = P(x,0), P2 = P(x,1), P3 = P(x,2), P4 = P(x,3);
    if (P1 < 0 || P2 < 0 || P3 < 0 || P4 < 0 ||
        P1 >= (int)nInd || P2 >= (int)nInd ||
        P3 >= (int)nInd || P4 >= (int)nInd) return;

    // Expected breeding values: one quarter of the four parental values
    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      const double eG = 0.25 * (GEBV(P1,ti) + GEBV(P2,ti) +
                               GEBV(P3,ti) + GEBV(P4,ti));
      results2(x, OFF_EG + ti) = eG;
      results2(x, OFF_SPV + ti) = eG;
    }

    // Select markers differing among the four parents. They are kept in
    // genetic-map order. A marker contributes nothing if all four agree.
    std::vector<std::vector<arma::uword>> chrDiff(chromosomeRanges.size());
    for (std::size_t c = 0; c < chromosomeRanges.size(); ++c) {
      const arma::uword startC = chromosomeRanges[c].first;
      const arma::uword endC   = chromosomeRanges[c].second;
      for (arma::uword mi = startC; mi <= endC; ++mi) {
        const double m1 = M_mat(P1,mi), m2 = M_mat(P2,mi);
        const double m3 = M_mat(P3,mi), m4 = M_mat(P4,mi);
        if (std::abs(m1-m2) > 1e-12 || std::abs(m3-m4) > 1e-12 ||
            std::abs(m1+m2-m3-m4) > 1e-12) {
          chrDiff[c].push_back(mi);
        }
      }
    }

    // Calculate only the upper triangle of the trait covariance matrix.
    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      for (arma::uword tj = ti; tj < numTrait; ++tj) {
        if (!covariance && ti != tj) continue;

        double SigmaSqP1P2 = 0.0;

        for (std::size_t c = 0; c < chromosomeRanges.size(); ++c) {
          const std::vector<arma::uword>& differing = chrDiff[c];
          if (differing.empty()) continue;

          // One accumulator per exponential component and contrast/trait.
          std::vector<double> h12_ti(K, 0.0), h12_tj(K, 0.0);
          std::vector<double> h34_ti(K, 0.0), h34_tj(K, 0.0);
          std::vector<double> hB_ti(K, 0.0),  hB_tj(K, 0.0);

          for (std::size_t a = 0; a < differing.size(); ++a) {
            const arma::uword mi = differing[a];

            // Weighted contrasts in the appendix's four-way equation.
            const double delta12 = M_mat(P1,mi) - M_mat(P2,mi);
            const double delta34 = M_mat(P3,mi) - M_mat(P4,mi);
            const double deltaB  = M_mat(P1,mi) + M_mat(P2,mi)
                                 - M_mat(P3,mi) - M_mat(P4,mi);

            const double v12_ti = delta12 * U_mat(mi,ti);
            const double v12_tj = delta12 * U_mat(mi,tj);
            const double v34_ti = delta34 * U_mat(mi,ti);
            const double v34_tj = delta34 * U_mat(mi,tj);
            const double vB_ti  = deltaB  * U_mat(mi,ti);
            const double vB_tj  = deltaB  * U_mat(mi,tj);

            // Correlation from the previous *retained* marker.
            double rho = 1.0;
            if (a > 0) {
              rho = std::exp(-2.0 *
                (markerPos[mi] - markerPos[differing[a-1]]));
            }

            double H12_ti = 0.0, H12_tj = 0.0;
            double H34_ti = 0.0, H34_tj = 0.0;
            double HB_ti = 0.0,  HB_tj = 0.0;
            double rhoPower = rho;

            for (int k = 0; k < K; ++k) {
              // Attenuate all earlier contributions by rho^(k+1).
              h12_ti[k] *= rhoPower; h12_tj[k] *= rhoPower;
              h34_ti[k] *= rhoPower; h34_tj[k] *= rhoPower;
              hB_ti[k]  *= rhoPower; hB_tj[k]  *= rhoPower;

              H12_ti += coeffW[k] * h12_ti[k];
              H12_tj += coeffW[k] * h12_tj[k];
              H34_ti += coeffW[k] * h34_ti[k];
              H34_tj += coeffW[k] * h34_tj[k];
              HB_ti  += coeffB[k] * hB_ti[k];
              HB_tj  += coeffB[k] * hB_tj[k];

              // Current marker enters the accumulator AFTER its contribution.
              h12_ti[k] += v12_ti; h12_tj[k] += v12_tj;
              h34_ti[k] += v34_ti; h34_tj[k] += v34_tj;
              hB_ti[k]  += vB_ti;  hB_tj[k]  += vB_tj;

              rhoPower *= rho;
            }

            // Diagonal marker term and all earlier marker-pair terms.
            // coeffW/B already contain the appendix's outer factor 1/4.
            SigmaSqP1P2 +=
                diagonalW * (v12_ti*v12_tj + v34_ti*v34_tj)
              + diagonalB * vB_ti*vB_tj
              + v12_ti*H12_tj + v12_tj*H12_ti
              + v34_ti*H34_tj + v34_tj*H34_ti
              + vB_ti*HB_tj   + vB_tj*HB_ti;
          }
        }

        const arma::uword kcol = tri_u_idx_incl(ti, tj);
        results1(x, kcol) = SigmaSqP1P2;
        if (ti == tj) {
          const double eG = results2(x, OFF_EG + ti);
          const double var = (SigmaSqP1P2 < 0.0 && SigmaSqP1P2 > -1e-12)
                           ? 0.0 : SigmaSqP1P2;
          results2(x, OFF_VAR + ti) = var;
          results2(x, OFF_SPV + ti) = eG + intensity * std::sqrt(var);
        }
      }
    }
  });

  // Reconstruct symmetric trait covariance matrices and optional index.
  if (covariance) {
    std::vector<arma::mat> covs(numCrosses);

    ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
      arma::uword x = static_cast<arma::uword>(xi);
      arma::mat V(numTrait, numTrait, arma::fill::zeros);
      for (arma::uword ti = 0; ti < numTrait; ++ti) {
        for (arma::uword tj = ti; tj < numTrait; ++tj) {
          const arma::uword k = tri_u_idx_incl(ti, tj);
          const double value = results1(x, k);
          V(ti,tj) = value;
          V(tj,ti) = value;
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
        const double sigma_gains = arma::as_scalar(weights_vec.t()*V*weights_vec);
        const double index = arma::as_scalar(
          results2.row(x).cols(0, numTrait-1) * weights_vec);
        const double spvi = index + intensity * std::sqrt(std::max(0.0, sigma_gains));
        results2(x,3*numTrait)     = index;
        results2(x,3*numTrait + 1) = sigma_gains;
        results2(x,3*numTrait + 2) = spvi;
      }
    });

    Rcpp::List out(numCrosses);
    for (arma::uword x = 0; x < numCrosses; ++x) out[x] = Rcpp::wrap(covs[x]);
    return Rcpp::List::create(
      Rcpp::Named("cross_values") = results2,
      Rcpp::Named("covariances") = out,
      Rcpp::Named("check_psd") = any_psd.load(std::memory_order_relaxed)
    );
  }
  return Rcpp::wrap(results2);
}
