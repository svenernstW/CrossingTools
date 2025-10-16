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
// [[Rcpp::export]]
SEXP cpp_calculate_desired_gains(const NumericMatrix& A,
                                     const NumericMatrix& V,
                                     const NumericMatrix& approxV,
                                     const NumericVector& gains,
                                     bool useMargV = false,
                                     bool useV = false,
                                     bool useapproxV = false,
                                     int nThreads = 4) {
#ifdef _OPENMP
  omp_set_dynamic(1);
  omp_set_num_threads(nThreads);
#endif
  const arma::uword numGeno   = A.nrow();
  const arma::uword numTrait     = A.ncol();
  arma::mat A_mat = as<arma::mat>(A);   // (numGeno × numTrait)
  arma::mat V_mat = as<arma::mat>(V);   // (numGeno × numTrait)
  arma::mat Vapprox_mat = as<arma::mat>(approxV);   // (numGeno × numTrait)
  
  arma::vec gains_vec = as<arma::vec>(gains);  // length == numTrait
  arma::vec results1(numGeno, arma::fill::zeros);
  
  if (useMargV) {
    arma::mat Vmarg(numTrait, numTrait, arma::fill::zeros);
#pragma omp parallel for schedule(dynamic)
    
    for (arma::uword gi = 0; gi < numGeno; ++gi) {
      arma::span idx = arma::span(gi * numTrait, gi * numTrait + numTrait - 1);
      
      Vmarg += V_mat(idx,idx);
      
    }
    
    Vmarg = Vmarg/numGeno;
   
    arma::mat L;
    bool spd = arma::chol(L, Vmarg);
    if (!spd) {
      arma::vec eval; arma::mat evec;
      arma::eig_sym(eval, evec, Vmarg);                 // O(n^3) but n is tiny
      double tol = std::max(1e-12, 1e-8 * eval.max());
      for (auto& l : eval) if (l < tol) l = tol;    // clamp
      Vmarg = evec * arma::diagmat(eval) * evec.t();    // PSD → SPD (since tol>0)
      arma::chol(L, Vmarg);                             // must succeed now
    }
    
    arma::vec w;
    arma::solve(w, Vmarg, gains_vec, arma::solve_opts::likely_sympd);

      
    results1 = A_mat * w;
  
    

    
    
  }
  
  if (useV) {
    arma::mat Vtemp(numTrait, numTrait, arma::fill::zeros);
    for (arma::uword gi = 0; gi < numGeno; ++gi) {
      arma::span idx = arma::span(gi * numTrait, gi * numTrait + numTrait - 1);
      Vtemp =  V_mat(idx,idx);
      
      arma::mat L;
      bool spd = arma::chol(L, Vtemp);
      if (!spd) {
        arma::vec eval; arma::mat evec;
        arma::eig_sym(eval, evec, Vtemp);                 // O(n^3) but n is tiny
        double tol = std::max(1e-12, 1e-8 * eval.max());
        for (auto& l : eval) if (l < tol) l = tol;    // clamp
        Vtemp = evec * arma::diagmat(eval) * evec.t();    // PSD → SPD (since tol>0)
        arma::chol(L, Vtemp);                             // must succeed now
      }
      
      arma::vec w;
      arma::solve(w, Vtemp, gains_vec, arma::solve_opts::likely_sympd);
      
      
      results1(gi) = dot(A_mat.row(gi), w);
      
    }
  }
  if (useapproxV) {
    
    arma::mat L;
    bool spd = arma::chol(L, Vapprox_mat);
    if (!spd) {
      arma::vec eval; arma::mat evec;
      arma::eig_sym(eval, evec, Vapprox_mat);                 // O(n^3) but n is tiny
      double tol = std::max(1e-12, 1e-8 * eval.max());
#pragma omp parallel for schedule(dynamic)
    for (auto& l : eval) if (l < tol) l = tol;    // clamp
      Vapprox_mat = evec * arma::diagmat(eval) * evec.t();    // PSD → SPD (since tol>0)
      arma::chol(L, Vapprox_mat);                             // must succeed now
    }
    
    arma::vec w;
    arma::solve(w, Vapprox_mat, gains_vec, arma::solve_opts::likely_sympd);
    
    
    results1 = A_mat * w;
    
    
    
    
    
  }
  return Rcpp::wrap(results1); 
  
}