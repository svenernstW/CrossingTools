// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <atomic>
#include <algorithm>
#include <limits>
#include <utility>

#include "parallel_backend.h"  // OpenMP when available, else RcppParallel

using namespace Rcpp;
using namespace arma;

// Haldane correlation of gametic origin between two markers.
// This is 1 - 2*r and can be propagated between retained markers.
inline double qjk_wolfe(double x, double y) {
  double diff = std::abs(x - y);
  double rcf  = 0.5 * (1.0 - std::exp(-2.0 * diff));
  return 1.0 - 2.0 * rcf;
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
                                    const NumericVector& p,
                                    bool covariance = false,
                                    bool calcindex = false,
                                    int nThreads = 4) {

  ct_set_threads(nThreads);
  std::atomic<bool> any_psd(false);

  const arma::uword numCrosses   = Crosses.nrow();
  const arma::uword numMarkers   = Hap1.ncol();
  const arma::uword numTrait     = U.ncol();
  const arma::uword numTraitComb = numTrait * (numTrait + 1) / 2;

  // Convert inputs to Armadillo
  arma::mat Hap1_mat = as<arma::mat>(Hap1);
  arma::mat Hap2_mat = as<arma::mat>(Hap2);
  arma::mat M_mat    = Hap1_mat + Hap2_mat;  // allele dosage 0..2
  arma::mat U_mat    = as<arma::mat>(U);
  arma::mat D_mat    = as<arma::mat>(D);
  arma::vec weights_vec = as<arma::vec>(weights);
  arma::rowvec p_vec    = as<arma::rowvec>(p);

  const arma::uword nInd = M_mat.n_rows;
  if (Hap2.nrow() != Hap1.nrow() || Hap2.ncol() != Hap1.ncol() ||
      U.nrow() != static_cast<int>(numMarkers) ||
      D.nrow() != static_cast<int>(numMarkers) ||
      D.ncol() != static_cast<int>(numTrait) ||
      p.size() != static_cast<R_xlen_t>(numMarkers) ||
      weights.size() != static_cast<R_xlen_t>(numTrait) ||
      Crosses.ncol() != 2 || numTrait == 0) {
    Rcpp::stop("Incompatible input dimensions");
  }
  if (numCrosses > static_cast<arma::uword>(std::numeric_limits<int>::max())) {
    Rcpp::stop("Too many crosses for the parallel backend");
  }

  // Genomic breeding values: (M - 2p) * U, without allocating W_mat
  arma::mat GEBV = M_mat * U_mat;
  arma::rowvec referenceMean = 2.0 * (p_vec * U_mat);
  GEBV.each_row() -= referenceMean;

  const arma::uword OFF_EG     = 0;
  const arma::uword OFF_ETG    = numTrait;
  const arma::uword OFF_VAR_A  = 2 * numTrait;
  const arma::uword OFF_SPV    = 3 * numTrait;
  const arma::uword OFF_VAR_D  = 4 * numTrait;
  const arma::uword OFF_VAR_T  = 5 * numTrait;
  const arma::uword OFF_TSPV   = 6 * numTrait;

  // Store ordered map positions and chromosome ranges; no quadratic RC matrix.
  arma::vec markerPos(numMarkers);
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
        Rcpp::stop("Map positions must be finite and ordered within chromosomes");
      }
      markerPos[startIdx + j] = pos;
    }
    chromosomeRanges.emplace_back(startIdx, startIdx + n - 1);
    startIdx += n;
  }
  if (startIdx != numMarkers) {
    Rcpp::stop("genMap and haplotype matrices have different marker counts");
  }

  // Copy parent indices (1-based in R to 0-based in C++)
  arma::Col<int> P1_idx(numCrosses), P2_idx(numCrosses);
  for (arma::uword x = 0; x < numCrosses; ++x) {
    P1_idx[x] = static_cast<int>(Crosses(x, 0)) - 1;
    P2_idx[x] = static_cast<int>(Crosses(x, 1)) - 1;
  }

  // Same result matrices and offsets as the original Wolfe function
  arma::mat results1A(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results1D(numCrosses, numTraitComb, arma::fill::zeros);
  // Cov(A,D) + Cov(D,A); on the diagonal this is 2*Cov(A,D)
  arma::mat results1AD(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results2(numCrosses, numTrait * 7 + 7, arma::fill::zeros);

  auto tri_u_idx_incl = [numTrait](arma::uword ti, arma::uword tj) -> arma::uword {
    return ti * numTrait - (ti * (ti - 1)) / 2 + (tj - ti);
  };

  // Parallel over crosses
  ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
    const arma::uword x = static_cast<arma::uword>(xi);
    const int P1 = P1_idx[x];
    const int P2 = P2_idx[x];
    if (P1 < 0 || P2 < 0 || P1 >= static_cast<int>(nInd) ||
        P2 >= static_cast<int>(nInd)) return;

    // Expected breeding values and statistical dominance deviations.
    // E[Z_k] = -0.5 * (dosage_P1 - 2p) * (dosage_P2 - 2p).
    arma::vec meanD(numMarkers);
    for (arma::uword mi = 0; mi < numMarkers; ++mi) {
      const double w1 = M_mat(P1, mi) - 2.0 * p_vec(mi);
      const double w2 = M_mat(P2, mi) - 2.0 * p_vec(mi);
      meanD[mi] = -0.5 * w1 * w2;
    }
    arma::rowvec eD = meanD.t() * D_mat;

    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      const double eG = 0.5 * (GEBV(P1, ti) + GEBV(P2, ti));
      results2(x, OFF_EG + ti)  = eG;
      results2(x, OFF_ETG + ti) = eG + eD[ti];
    }

    // Process each chromosome independently. Skip any locus at which
    // both parents are homozygous: neither gamete segregates there.
    for (const auto& chr : chromosomeRanges) {
      const arma::uword startC = chr.first;
      const arma::uword endC   = chr.second;

      std::vector<arma::uword> chrDiff;
      for (arma::uword mi = startC; mi <= endC; ++mi) {
        if (Hap1_mat(P1, mi) != Hap2_mat(P1, mi) ||
            Hap1_mat(P2, mi) != Hap2_mat(P2, mi)) {
          chrDiff.push_back(mi);
        }
      }
      if (chrDiff.empty()) continue;

      // Five running components (one per trait):
      // additive gamete 1, additive gamete 2,
      // dominance gamete 1, dominance gamete 2, dominance interaction.
      // First four use exp(-2d); interaction uses exp(-4d).
      arma::rowvec hA1(numTrait, arma::fill::zeros);
      arma::rowvec hA2(numTrait, arma::fill::zeros);
      arma::rowvec hD1(numTrait, arma::fill::zeros);
      arma::rowvec hD2(numTrait, arma::fill::zeros);
      arma::rowvec hD12(numTrait, arma::fill::zeros);

      arma::rowvec a1(numTrait), a2(numTrait);
      arma::rowvec d1(numTrait), d2(numTrait), d12(numTrait);

      for (std::size_t ii = 0; ii < chrDiff.size(); ++ii) {
        const arma::uword mi = chrDiff[ii];

        // Propagate previous segregating markers to the current position.
        if (ii > 0) {
          const double rho = qjk_wolfe(markerPos[mi],
                                       markerPos[chrDiff[ii - 1]]);
          hA1 *= rho;
          hA2 *= rho;
          hD1 *= rho;
          hD2 *= rho;
          hD12 *= rho * rho;
        }

        // Haplotypic contrasts and expected transmitted alleles
        const double delta1 = Hap1_mat(P1, mi) - Hap2_mat(P1, mi);
        const double delta2 = Hap1_mat(P2, mi) - Hap2_mat(P2, mi);
        const double EA = 0.5 * M_mat(P1, mi);
        const double EB = 0.5 * M_mat(P2, mi);

        // Z - E[Z] = 2(p-EB)(A-EA) + 2(p-EA)(B-EB)
        //            - 2(A-EA)(B-EB).
        const double deltaD1 = 2.0 * (p_vec(mi) - EB) * delta1;
        const double deltaD2 = 2.0 * (p_vec(mi) - EA) * delta2;
        const double deltaD12 = delta1 * delta2;

        // Marker effects for each of the five covariance components.
        for (arma::uword ti = 0; ti < numTrait; ++ti) {
          a1[ti] = delta1 * U_mat(mi, ti);
          a2[ti] = delta2 * U_mat(mi, ti);
          d1[ti] = deltaD1 * D_mat(mi, ti);
          d2[ti] = deltaD2 * D_mat(mi, ti);
          d12[ti] = deltaD12 * D_mat(mi, ti);
        }

        // Only the upper triangle is evaluated. Off-diagonal trait
        // covariances are skipped entirely when covariance == false.
        for (arma::uword ti = 0; ti < numTrait; ++ti) {
          const arma::uword endTrait = covariance ? numTrait : ti + 1;
          for (arma::uword tj = ti; tj < endTrait; ++tj) {
            const arma::uword kcol = tri_u_idx_incl(ti, tj);

            // Additive-additive: two independent gametic contributions.
            results1A(x, kcol) += 0.25 * (
                a1[ti] * a1[tj] + a2[ti] * a2[tj]
              + a1[ti] * hA1[tj] + a1[tj] * hA1[ti]
              + a2[ti] * hA2[tj] + a2[tj] * hA2[ti]
            );

            // Dominance-dominance: two exp(-2d) terms and one exp(-4d) term.
            results1D(x, kcol) += 0.25 * (
                d1[ti] * d1[tj] + d2[ti] * d2[tj]
              + d12[ti] * d12[tj]
              + d1[ti] * hD1[tj] + d1[tj] * hD1[ti]
              + d2[ti] * hD2[tj] + d2[tj] * hD2[ti]
              + d12[ti] * hD12[tj] + d12[tj] * hD12[ti]
            );

            // Cov(additive, dominance) + Cov(dominance, additive).
            results1AD(x, kcol) += 0.25 * (
                a1[ti] * d1[tj] + d1[ti] * a1[tj]
              + a2[ti] * d2[tj] + d2[ti] * a2[tj]
              + a1[ti] * hD1[tj] + d1[tj] * hA1[ti]
              + d1[ti] * hA1[tj] + a1[tj] * hD1[ti]
              + a2[ti] * hD2[tj] + d2[tj] * hA2[ti]
              + d2[ti] * hA2[tj] + a2[tj] * hD2[ti]
            );
          }
        }

        // The current marker becomes part of the previous-marker sums.
        hA1 += a1;
        hA2 += a2;
        hD1 += d1;
        hD2 += d2;
        hD12 += d12;
      }
    }

    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      const arma::uword kdiag = tri_u_idx_incl(ti, ti);
      const double varA = results1A(x, kdiag);
      const double varD = results1D(x, kdiag);
      const double varAD = results1AD(x, kdiag);
      const double varT = std::max(0.0, varA + varD + varAD);
      const double eG = results2(x, OFF_EG + ti);
      const double eTG = results2(x, OFF_ETG + ti);

      results2(x, OFF_VAR_A + ti) = varA;
      results2(x, OFF_VAR_D + ti) = varD;
      results2(x, OFF_VAR_T + ti) = varT;
      results2(x, OFF_SPV + ti) = eG + intensity * std::sqrt(std::max(0.0, varA));
      results2(x, OFF_TSPV + ti) = eTG + intensity * std::sqrt(varT);
    }
  });

  // Same output structure and index calculations as the original implementation.
  if (covariance) {
    std::vector<arma::mat> covsA(numCrosses);
    std::vector<arma::mat> covsD(numCrosses);
    std::vector<arma::mat> covsAD(numCrosses);

    const arma::uword OFF_IDX_A     = 7 * numTrait + 0;
    const arma::uword OFF_IDX_T     = 7 * numTrait + 1;
    const arma::uword OFF_VARIDX_A  = 7 * numTrait + 2;
    const arma::uword OFF_SPVI_A    = 7 * numTrait + 3;
    const arma::uword OFF_VARIDX_D  = 7 * numTrait + 4;
    const arma::uword OFF_VARIDX_T  = 7 * numTrait + 5;
    const arma::uword OFF_SPVI_T    = 7 * numTrait + 6;

    ct_parallel_for(0, static_cast<int>(numCrosses), [&](int xi) {
      arma::uword x = static_cast<arma::uword>(xi);
      arma::mat GA(numTrait, numTrait, arma::fill::zeros);
      arma::mat GD(numTrait, numTrait, arma::fill::zeros);
      arma::mat GAD(numTrait, numTrait, arma::fill::zeros);

      for (arma::uword ti = 0; ti < numTrait; ++ti) {
        for (arma::uword tj = ti; tj < numTrait; ++tj) {
          const arma::uword k = tri_u_idx_incl(ti, tj);
          GA(ti, tj) = GA(tj, ti) = results1A(x, k);
          GD(ti, tj) = GD(tj, ti) = results1D(x, k);
          GAD(ti, tj) = GAD(tj, ti) = results1AD(x, k);
        }
      }
      covsA[x] = GA;
      covsD[x] = GD;
      covsAD[x] = GAD;

      if (calcindex) {
        arma::mat VA = GA;
        arma::mat VD = GD;
        arma::mat VT = GA + GD + GAD;
        arma::mat L;

        if (!arma::chol(L, VA)) {
          any_psd.store(true, std::memory_order_relaxed);
          arma::vec eval;
          arma::mat evec;
          arma::eig_sym(eval, evec, VA);
          eval.transform([](double v) { return (v < 0.0) ? 0.0 : v; });
          VA = evec * arma::diagmat(eval) * evec.t();
        }
        if (!arma::chol(L, VD)) {
          any_psd.store(true, std::memory_order_relaxed);
          arma::vec eval;
          arma::mat evec;
          arma::eig_sym(eval, evec, VD);
          eval.transform([](double v) { return (v < 0.0) ? 0.0 : v; });
          VD = evec * arma::diagmat(eval) * evec.t();
        }
        if (!arma::chol(L, VT)) {
          any_psd.store(true, std::memory_order_relaxed);
          arma::vec eval;
          arma::mat evec;
          arma::eig_sym(eval, evec, VT);
          eval.transform([](double v) { return (v < 0.0) ? 0.0 : v; });
          VT = evec * arma::diagmat(eval) * evec.t();
        }

        const double index_A = arma::as_scalar(
          results2.row(x).cols(OFF_EG, OFF_EG + numTrait - 1) * weights_vec);
        const double index_T = arma::as_scalar(
          results2.row(x).cols(OFF_ETG, OFF_ETG + numTrait - 1) * weights_vec);
        const double var_index_A = std::max(0.0,
          arma::as_scalar(weights_vec.t() * VA * weights_vec));
        const double var_index_D = std::max(0.0,
          arma::as_scalar(weights_vec.t() * VD * weights_vec));
        const double var_index_T = std::max(0.0,
          arma::as_scalar(weights_vec.t() * VT * weights_vec));

        results2(x, OFF_IDX_A) = index_A;
        results2(x, OFF_IDX_T) = index_T;
        results2(x, OFF_VARIDX_A) = var_index_A;
        results2(x, OFF_SPVI_A) = index_A + intensity * std::sqrt(var_index_A);
        results2(x, OFF_VARIDX_D) = var_index_D;
        results2(x, OFF_VARIDX_T) = var_index_T;
        results2(x, OFF_SPVI_T) = index_T + intensity * std::sqrt(var_index_T);
      }
    });

    Rcpp::List out_A(numCrosses), out_D(numCrosses), out_AD(numCrosses);
    for (arma::uword x = 0; x < numCrosses; ++x) {
      out_A[x] = Rcpp::wrap(covsA[x]);
      out_D[x] = Rcpp::wrap(covsD[x]);
      out_AD[x] = Rcpp::wrap(covsAD[x]);
    }
    return Rcpp::List::create(
      Rcpp::Named("cross_values") = results2,
      Rcpp::Named("covA") = out_A,
      Rcpp::Named("covD") = out_D,
      Rcpp::Named("covAD") = out_AD,
      Rcpp::Named("check_psd") = any_psd.load(std::memory_order_relaxed)
    );
  }

  return Rcpp::wrap(results2);
}
