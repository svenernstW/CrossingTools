// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include "parallel_backend.h"  // OpenMP when available, else RcppParallel

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP cpp_calculate_desired_gains(const NumericMatrix& A,
                                 const NumericMatrix& V,
                                 const NumericMatrix& approxV,
                                 const NumericVector& gains,
                                 bool useMargV = false,
                                 bool useV = false,
                                 bool useapproxV = false,
                                 int nThreads = 4) {
  // portable thread setup
  ct_set_threads(nThreads);

  const arma::uword numGeno  = A.nrow();
  const arma::uword numTrait = A.ncol();

  arma::mat A_mat       = as<arma::mat>(A);        // (numGeno × numTrait)
  arma::mat V_mat       = as<arma::mat>(V);        // (numGeno*numTrait × numGeno*numTrait) block-diag by geno
  arma::mat Vapprox_mat = as<arma::mat>(approxV);  // (numTrait × numTrait)
  arma::vec gains_vec   = as<arma::vec>(gains);    // length == numTrait

  arma::vec results1(numGeno, arma::fill::zeros);

  // 1) Use average marginal V over genotypes
  if (useMargV) {
    arma::mat Vmarg(numTrait, numTrait, arma::fill::zeros);

    // Sum blocks serially (safe & numTrait is small; parallel reduction would add complexity)
    for (arma::uword gi = 0; gi < numGeno; ++gi) {
      arma::span idx = arma::span(gi * numTrait, gi * numTrait + numTrait - 1);
      Vmarg += V_mat(idx, idx);
    }
    Vmarg /= static_cast<double>(numGeno);

    arma::mat L;
    if (!arma::chol(L, Vmarg)) {
      arma::vec eval; arma::mat evec;
      arma::eig_sym(eval, evec, Vmarg);
      double tol = std::max(1e-12, 1e-8 * eval.max());
      for (auto &l : eval) if (l < tol) l = tol;
      Vmarg = evec * arma::diagmat(eval) * evec.t();
      arma::chol(L, Vmarg);
    }

    arma::vec w;
    arma::solve(w, Vmarg, gains_vec, arma::solve_opts::likely_sympd);
    results1 = A_mat * w;  // per-genotype index
  }

  // 2) Use each genotype's own V block
  if (useV) {
    // Each gi independent → parallelize safely
    ct_parallel_for(0, static_cast<int>(numGeno), [&](int gi_i) {
      arma::uword gi = static_cast<arma::uword>(gi_i);
      arma::span idx = arma::span(gi * numTrait, gi * numTrait + numTrait - 1);
      arma::mat Vtemp = V_mat(idx, idx);

      arma::mat L;
      if (!arma::chol(L, Vtemp)) {
        arma::vec eval; arma::mat evec;
        arma::eig_sym(eval, evec, Vtemp);
        double tol = std::max(1e-12, 1e-8 * eval.max());
        for (auto &l : eval) if (l < tol) l = tol;
        Vtemp = evec * arma::diagmat(eval) * evec.t();
        arma::chol(L, Vtemp);
      }

      arma::vec w;
      arma::solve(w, Vtemp, gains_vec, arma::solve_opts::likely_sympd);
      results1(gi) = arma::dot(A_mat.row(gi), w);
    });
  }

  // 3) Use a precomputed approximate V (same for all genotypes)
  if (useapproxV) {
    arma::mat L;
    if (!arma::chol(L, Vapprox_mat)) {
      arma::vec eval; arma::mat evec;
      arma::eig_sym(eval, evec, Vapprox_mat);
      double tol = std::max(1e-12, 1e-8 * eval.max());
      for (auto &l : eval) if (l < tol) l = tol;
      Vapprox_mat = evec * arma::diagmat(eval) * evec.t();
      arma::chol(L, Vapprox_mat);
    }

    arma::vec w;
    arma::solve(w, Vapprox_mat, gains_vec, arma::solve_opts::likely_sympd);
    results1 = A_mat * w;
  }

  return Rcpp::wrap(results1);
}
