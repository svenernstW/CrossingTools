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
  ct_set_threads(nThreads);

  const arma::uword numGeno  = A.nrow();
  const arma::uword numTrait = A.ncol();

  arma::mat A_mat       = as<arma::mat>(A);
  arma::mat V_mat       = as<arma::mat>(V);
  arma::mat Vapprox_mat = as<arma::mat>(approxV);
  arma::vec gains_vec   = as<arma::vec>(gains);

  arma::vec results1(numGeno, arma::fill::zeros);
  arma::vec w_out;          // only valid for useMargV/useapproxV
  bool have_w = false;

  if (useMargV) {
    arma::mat Vmarg(numTrait, numTrait, arma::fill::zeros);
    for (arma::uword gi = 0; gi < numGeno; ++gi) {
      arma::span idx(gi * numTrait, gi * numTrait + numTrait - 1);
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

    arma::solve(w_out, Vmarg, gains_vec, arma::solve_opts::likely_sympd);
    results1 = A_mat * w_out;
    have_w = true;
  }

  if (useV) {
    ct_parallel_for(0, static_cast<int>(numGeno), [&](int gi_i) {
      arma::uword gi = static_cast<arma::uword>(gi_i);
      arma::span idx(gi * numTrait, gi * numTrait + numTrait - 1);
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

      arma::vec w;  // thread-local (no race)
      arma::solve(w, Vtemp, gains_vec, arma::solve_opts::likely_sympd);
      results1(gi) = arma::dot(A_mat.row(gi), w);
    });
    // have_w stays false: no single global weight
  }

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

    arma::solve(w_out, Vapprox_mat, gains_vec, arma::solve_opts::likely_sympd);
    results1 = A_mat * w_out;
    have_w = true;
  }

  return Rcpp::List::create(
    Named("index")  = results1,
    Named("weight") = have_w ? Rcpp::wrap(w_out) : R_NilValue
  );
}
