// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include "parallel_backend.h"  // portable threading: OpenMP if present, else RcppParallel

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP cpp_u_from_from_g(const NumericMatrix& M,
                       const NumericMatrix& G,
                       const NumericMatrix& g,
                       double scalingFactor,
                       double tol = 1e-10,
                       bool LDvar = false,
                       const Nullable<DataFrame>& Grouping = R_NilValue,
                       bool calcPriorVcov = false,
                       const Nullable<NumericMatrix> PriorVcov = R_NilValue,
                       const Nullable<NumericMatrix> sigmasq = R_NilValue,
                       bool calcPosteriorVcov = false,
                       const Nullable<NumericMatrix> PEV = R_NilValue,
                       int nThreads = 4) {

  // portable thread setup
  ct_set_threads(nThreads);

  // Convert inputs to Armadillo
  mat M_mat = as<mat>(M);        // n x p
  mat G_mat = as<mat>(G);        // n x n
  mat g_mat = as<mat>(g);        // n x k
  int n = M_mat.n_rows;
  int p = M_mat.n_cols;
  int k = g_mat.n_cols;

  // Grouping (for LDvar)
  IntegerVector groups_R;
  if (LDvar) {
    if (Grouping.isNotNull()) {
      DataFrame df(Grouping.get());
      if (!df.containsElementNamed("chr"))
        stop("Grouping must contain a 'chr' column.");
      groups_R = df["chr"];
      if ((int)groups_R.size() != p)
        stop("Length of Grouping$chr must equal number of markers (ncol(M)).");
    } else {
      groups_R = IntegerVector(p, 1);
    }
  }

  // Check/regularize G, then invert
  mat L,U,Pmat;
  lu(L, U, Pmat, G_mat);
  if (U.is_empty() || U.diag().min() < tol) {
    Rcpp::Rcout << "G matrix near-singular; adding 1e-4 to diagonal.\n";
    G_mat.diag() += 1e-4;
  }
  mat Gi;
  bool ok = inv_sympd(Gi, G_mat);
  if (!ok) Gi = inv(G_mat);

  // Backsolve marker effects: mu = (scalingFactor * M' * Gi) * g
  mat C  = (scalingFactor * M_mat.t()) * Gi;   // p x n
  mat mu = C * g_mat;                          // p x k

  // -------- Prior Vcov (trait–marker) --------
  mat priorVcov_pk;  // (p*k) x (p*k); left empty if !calcPriorVcov
  if (calcPriorVcov) {
    if (sigmasq.isNull()) {
      stop("calcPriorVcov = TRUE requires 'sigmasq' (k x k) or scalar/vector.");
    }

    // Parse/normalize Sigma (k x k)
    mat s = as<mat>(NumericMatrix(sigmasq));
    mat Sigma;
    if ((int)s.n_rows == k && (int)s.n_cols == k) {
      Sigma = s;  // full trait covariance provided
    } else if ((int)s.n_elem == 1) {
      Sigma = mat(k, k, fill::eye) * s(0,0);          // scalar -> scalar * I_k
    } else if ((int)s.n_elem == k && (s.n_rows == 1 || s.n_cols == 1)) {
      Sigma = diagmat(vectorise(s));                  // length-k vector -> diag
    } else {
      stop("sigmasq must be (k x k), scalar, or length-k vector.");
    }

    if (PriorVcov.isNotNull()) {
      mat Prior = as<mat>(NumericMatrix(PriorVcov));
      if ((int)Prior.n_rows == p*k && (int)Prior.n_cols == p*k) {
        // Full prior provided; assume already scaled appropriately.
        priorVcov_pk = Prior;
      } else if ((int)Prior.n_rows == p && (int)Prior.n_cols == p) {
        // Marker prior base -> expand with traits
        // Match your earlier convention: divide by p
        priorVcov_pk = kron(Prior, Sigma) / (double)p;
      } else {
        stop("PriorVcov must be either p x p or (p*k) x (p*k).");
      }
    } else {
      // No marker prior provided -> identity on markers
      Rcpp::Rcout << "PriorVcov not supplied; using kron(I_p, sigmasq) / p.\n";
      priorVcov_pk = kron(eye<mat>(p, p), Sigma) / (double)p;
    }
  }

  // -------- Posterior Vcov: PEV is (n*k) x (n*k) --------
  // Build a (p*k) x (p*k) block matrix: each (t,u) block is C * PEV_tu * C'
  mat posteriorVcov; // will be (p*k) x (p*k)
  if (calcPosteriorVcov) {
    if (PEV.isNull()) {
      stop("calcPosteriorVcov = TRUE but PEV is NULL. Expecting (n*k) x (n*k) PEV from mixed model.");
    }
    mat PEV_big = as<mat>(NumericMatrix(PEV));
    if ((int)PEV_big.n_rows != n*k || (int)PEV_big.n_cols != n*k) {
      stop("PEV must be of size (n*traits) x (n*traits).");
    }

    // Allocate result
    posteriorVcov.set_size(p*k, p*k);
    posteriorVcov.zeros();

    // Parallelize over all k*k blocks safely (each block writes to disjoint submat)
    ct_parallel_for(0, k*k, [&](int idx) {
      int t = idx / k;
      int u = idx % k;
      int rt  = t * n;
      int cu  = u * n;
      int rmu = t * p;
      int cmu = u * p;
      mat B = C * PEV_big.submat(rt, cu, rt + n - 1, cu + n - 1) * C.t(); // p x p
      posteriorVcov.submat(rmu, cmu, rmu + p - 1, cmu + p - 1) = B;
    });
  } else {
    posteriorVcov.set_size(0,0);
  }

  // -------- LD variance by groups (unchanged logic) --------
  if (!LDvar) {
    return List::create(
      Named("mu_matrix")      = mu,                                                     // p x k
      Named("variance_LD")    = R_NilValue,
      Named("prior_Vcov")     = (priorVcov_pk.is_empty() ? R_NilValue : wrap(priorVcov_pk)),
      Named("posterior_Vcov") = (posteriorVcov.n_elem ? posteriorVcov : mat(0,0))      // (p*k) x (p*k)
    );
  }

  // Build unique groups and indices
  std::unordered_set<int> uniq_set(groups_R.begin(), groups_R.end());
  std::vector<int> uniq_groups(uniq_set.begin(), uniq_set.end());
  std::sort(uniq_groups.begin(), uniq_groups.end());
  const size_t n_groups = uniq_groups.size();

  std::vector<uvec> group_indices(n_groups);
  for (size_t gi = 0; gi < n_groups; ++gi) {
    int gval = uniq_groups[gi];
    std::vector<uword> idx; idx.reserve(p);
    for (int col = 0; col < p; ++col)
      if (groups_R[col] == gval) idx.push_back((uword)col);
      group_indices[gi] = uvec(idx);
  }

  std::vector<mat> M_subs(n_groups);
  std::vector<mat> mu_subs(n_groups);
  for (size_t gi = 0; gi < n_groups; ++gi) {
    const uvec& idx = group_indices[gi];
    M_subs[gi]  = M_mat.cols(idx); // n x p_i
    mu_subs[gi] = mu.rows(idx);    // p_i x k
  }

  List ld_list(k);
  for (int t = 0; t < k; ++t) {
    mat Gcov(n_groups, n_groups, fill::zeros);
    for (size_t i = 0; i < n_groups; ++i) {
      vec mui = mu_subs[i].col(t);
      mat Zi  = M_subs[i];
      for (size_t j = i; j < n_groups; ++j) {
        vec muj = mu_subs[j].col(t);
        mat Zj  = M_subs[j];
        double val = as_scalar( mui.t() * (Zi.t() * (Zj * muj)) );
        if (i == j) Gcov(i,i) = val / (double)n;
        else {
          double off = (2.0 * val) / (double)n;
          Gcov(i,j) = off; Gcov(j,i) = off;
        }
      }
    }
    ld_list[t] = Gcov;
  }
  CharacterVector ld_names(k);
  for (int t = 0; t < k; ++t) ld_names[t] = "trait_" + std::to_string(t+1);
  ld_list.attr("names") = ld_names;

  return List::create(
    Named("mu_matrix")      = mu,                                                     // p x k
    Named("variance_LD")    = (k==1 ? ld_list[0] : ld_list),                          // matrix if k==1, else list
    Named("prior_Vcov")     = (priorVcov_pk.is_empty() ? R_NilValue : wrap(priorVcov_pk)),
    Named("posterior_Vcov") = (posteriorVcov.n_elem ? posteriorVcov : mat(0,0)),      // (p*k) x (p*k)
    Named("group_levels")   = wrap(uniq_groups)
  );
}
