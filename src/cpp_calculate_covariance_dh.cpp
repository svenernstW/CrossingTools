// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <atomic>
#include <limits>
#include "parallel_backend.h"   // OpenMP when available, else RcppParallel

using namespace Rcpp;
using namespace arma;

// Explicit DH recombination kernel (reference formula).
// Breeding order: t random-mating rounds, then s selfing rounds, then DH.
// The variance calculation below uses its exponential expansion instead.
inline double qjk_dh(double x, double y, int t, int s) {
  double diff = std::abs(x - y);
  double rcf  = 0.5 * (1.0 - std::exp(-2.0 * diff));

  // The first random-mating round of a uniform F1 gives the F2 distribution.
  // It is therefore counted in the DH-generation index below.
  int randomRounds = (t > 0) ? t - 1 : 0;
  int dhGeneration = s + 1 + ((t > 0) ? 1 : 0);

  // Kernel 1: remaining rounds of random intermating
  double randomMatingKernel = std::pow(1.0 - rcf, randomRounds);

  // Kernel 2: selfing followed by DH production (Lehermeier DH model)
  double selfingDhKernel =
    1.0 - (4.0 * rcf / (1.0 + 2.0 * rcf)) *
    (1.0 - std::pow(0.5 * (1.0 - 2.0 * rcf), dhGeneration));

  return 0.5 + 0.5 * randomMatingKernel * selfingDhKernel;
}

// [[Rcpp::export]]
SEXP cpp_calculate_covariance_dh(const NumericMatrix& Crosses,
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
  // portable thread setup
  ct_set_threads(nThreads);
  std::atomic<bool> any_psd(false);

  const arma::uword numCrosses   = Crosses.nrow();
  const arma::uword numMarkers   = M.ncol();
  const arma::uword numTrait     = U.ncol();
  const arma::uword numTraitComb = numTrait * (numTrait + 1) / 2;

  // Convert inputs to Armadillo
  arma::mat M_mat = as<arma::mat>(M);   // genotype dosage 0..2
  arma::mat U_mat = as<arma::mat>(U);   // average substitution effects alpha

  arma::vec weights_vec =
    as<arma::vec>(weights);

  // Allele frequencies in the reference population
  arma::rowvec p_vec = as<arma::rowvec>(p);

  // Centred additive genotype code W = M - 2p
  arma::mat W_mat = M_mat;
  W_mat.each_row() -= 2.0 * p_vec;

  // Genomic breeding values
  arma::mat GEBV =
    W_mat * U_mat;

  const arma::uword OFF_EG  = 0;
  const arma::uword OFF_VAR = numTrait;
  const arma::uword OFF_SPV = 2 * numTrait;

  // Precompute exponential kernel weights instead of the full QJK matrix.
  // The random-mating and selfing/DH kernels above are expanded using
  // z = 1 - 2*r = exp(-2*distance):
  //   randomMatingKernel = ((1+z)/2)^randomRounds
  //   selfingDhKernel    = sum_{k=1}^{g-1} z^k/2^k + 2*z^g/2^g
  // with g = dhGeneration. Their product has K = t+s+1 components.
  if (t < 0 || s < 0) {
    Rcpp::stop("t and s must be non-negative integers");
  }
  if (t > std::numeric_limits<int>::max() - s - 1) {
    Rcpp::stop("Too many generations");
  }

  const int randomRounds = (t > 0) ? t - 1 : 0;
  const int dhGeneration = s + 1 + ((t > 0) ? 1 : 0);
  const int K = randomRounds + dhGeneration;

  // Binomial weights of the random-mating kernel:
  // ((1 + z) / 2)^randomRounds
  std::vector<double> randomWeights(randomRounds + 1);
  randomWeights[0] = std::pow(0.5, randomRounds);
  for (int j = 1; j <= randomRounds; ++j) {
    randomWeights[j] = randomWeights[j - 1] *
      (randomRounds - j + 1) / j;
  }

  // Multiply the random-mating kernel by the selfing/DH kernel:
  // sum_{k=1}^{g-1} z^k/2^k + 2*z^g/2^g, where g = dhGeneration.
  // coeff[k-1] multiplies exp(-2*k*distance), including the factor 0.25.
  std::vector<double> coeff(K, 0.0);
  for (int j = 0; j <= randomRounds; ++j) {
    for (int k = 1; k <= dhGeneration; ++k) {
      double selfingWeight = std::pow(0.5, k);
      if (k == dhGeneration) selfingWeight *= 2.0;
      coeff[j + k - 1] += 0.25 * randomWeights[j] * selfingWeight;
    }
  }

  // Precompute genetic positions for each chromosome (no QJK matrix needed)
  arma::vec markerPos(numMarkers);
  arma::uword startIdx = 0;
  std::vector<std::pair<arma::uword, arma::uword>> chromosomeRanges;

  for (R_xlen_t c = 0; c < genMap.size(); ++c) {
    NumericMatrix gm = as<NumericMatrix>(genMap[c]);  // col 0 = position in Morgans
    const arma::uword n = gm.nrow();

    if (n == 0) continue;
    if (gm.ncol() < 1 || startIdx + n > numMarkers) {
      Rcpp::stop("Invalid genetic map dimensions");
    }
    for (arma::uword j = 0; j < n; ++j) {
      if (!std::isfinite(gm(j, 0)) ||
          (j > 0 && gm(j, 0) < gm(j - 1, 0))) {
        Rcpp::stop("Marker positions must be finite and ordered within chromosomes");
      }
      markerPos[startIdx + j] = gm(j, 0);
    }
    chromosomeRanges.emplace_back(startIdx, startIdx + n - 1);
    startIdx += n;
  }
  if (startIdx != numMarkers) {
    Rcpp::stop("genMap and M contain different marker counts");
  }

  // Number of parents (marker effects are applied only where needed)
  const arma::uword nInd = M_mat.n_rows;

  // Copy parent indices out of Crosses (1-based in R -> 0-based here)
  arma::Col<int> P1_idx(numCrosses), P2_idx(numCrosses);
  for (arma::uword x = 0; x < numCrosses; ++x) {
    P1_idx[x] = static_cast<int>(Crosses(x, 0)) - 1;
    P2_idx[x] = static_cast<int>(Crosses(x, 1)) - 1;
  }

  // Results
  arma::mat results1(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results2(numCrosses, numTrait * 3 + 3, arma::fill::zeros);

  // Helper to map (ti, tj) with 0 <= ti <= tj < numTrait to column index
  auto tri_u_idx_incl = [numTrait](arma::uword ti, arma::uword tj) -> arma::uword {
    return ti * numTrait - (ti * (ti - 1)) / 2 + (tj - ti);
  };

  // Parallel over crosses
  ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
    R_xlen_t x = static_cast<R_xlen_t>(xi);
    const int P1 = P1_idx[x];
    const int P2 = P2_idx[x];
    if (P1 < 0 || P2 < 0 || P1 >= (int)nInd || P2 >= (int)nInd) return;

    // Markers that differ between the two parents (exact compare; use tolerance if needed)
    arma::uvec differing = find(abs(M_mat.row(P1) - M_mat.row(P2)) > 1e-12);
    if (differing.n_elem == 0) {
      for (arma::uword ti = 0; ti < numTrait; ++ti) {
        double G1 = GEBV(P1, ti);
        double G2 = GEBV(P2, ti);
        double eG = 0.5 * (G1 + G2);

        results2(x, OFF_EG  + ti) = eG;
        results2(x, OFF_VAR + ti) = 0;
        results2(x, OFF_SPV + ti) = eG;
      }
      return;
    }
    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      for (arma::uword tj = ti; tj < numTrait; ++tj) {
        if (!covariance && ti != tj) continue;

        double SigmaSqP1P2 = 0.0;  // reset per trait pair

        // Accumulate by chromosome
        for (const auto& chr : chromosomeRanges) {
          const arma::uword startC = chr.first;
          const arma::uword endC   = chr.second;

          // differing markers within [startC, endC]
          arma::uvec idx    = find( (differing >= startC) % (differing <= endC) );
          if (idx.n_elem == 0) continue;
          arma::uvec chrDiff = differing.elem(idx);

          // One running accumulator per exponential kernel and trait.
          // These are reset for each chromosome and trait pair.
          arma::vec h_ti(K, arma::fill::zeros);
          arma::vec h_tj(K, arma::fill::zeros);

          const arma::uword mlen = chrDiff.n_elem;
          for (arma::uword a = 0; a < mlen; ++a) {
            const arma::uword mi = chrDiff[a];
            // Calculate marker contributions directly, without the P_mu cube.
            const double delta = M_mat(P1, mi) - M_mat(P2, mi);
            const double di_ti = delta * U_mat(mi, ti);
            const double di_tj = delta * U_mat(mi, tj);

            // Correlation between consecutive DIFFERING markers.
            // Skipped non-segregating markers are accounted for by the distance.
            double rho = 1.0;
            if (a > 0) {
              rho = std::exp(-2.0 *
                (markerPos[mi] - markerPos[chrDiff[a - 1]]));
            }

            // Replace the inner marker-pair loop by K linear updates.
            double H_ti = 0.0;
            double H_tj = 0.0;
            double rhoPower = rho;
            for (int k = 0; k < K; ++k) {
              h_ti[k] *= rhoPower;
              h_tj[k] *= rhoPower;

              H_ti += coeff[k] * h_ti[k];
              H_tj += coeff[k] * h_tj[k];

              h_ti[k] += di_ti;
              h_tj[k] += di_tj;
              rhoPower *= rho;
            }

            // Diagonal marker term plus both cross-marker/trait contributions
            SigmaSqP1P2 += 0.25 * di_ti * di_tj
              + di_ti * H_tj + di_tj * H_ti;
          }
        }

        const arma::uword k = tri_u_idx_incl(ti, tj);
        results1(x, k) = SigmaSqP1P2;

        if (ti == tj) {
          double G1 = GEBV(P1, ti);
          double G2 = GEBV(P2, ti);
          double eG = 0.5 * (G1 + G2);

          results2(x, OFF_EG  + ti) = eG;
          results2(x, OFF_VAR + ti) = SigmaSqP1P2;
          results2(x, OFF_SPV + ti) = eG + intensity * std::sqrt(SigmaSqP1P2);
        }
      }
    }
  });

  // If requested, convert each row to an nTrait x nTrait symmetric matrix
  if (covariance) {
    std::vector<arma::mat> covs(numCrosses);

    ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
      arma::uword x = static_cast<arma::uword>(xi);
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

      if (calcindex) {
        arma::mat V = covs[x];
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

    // wrap covariances into an R list
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
