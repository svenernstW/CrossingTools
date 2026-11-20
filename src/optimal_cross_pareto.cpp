// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <numeric>
#include <limits>
#include <sstream>
#include <cstdint>

#include "parallel_backend.h"  // OpenMP when available, else RcppParallel

using namespace Rcpp;

// ------------------------------
// Helpers
// ------------------------------
// Parental contributions for a crossing plan (indices are 0-based)
static arma::vec calcContr(const arma::umat& crosses, const arma::uword& nInd,
                           const arma::uword& nCross) {

  double val = 1.0 / (2.0 * static_cast<double>(nCross));
  arma::vec x(nInd, arma::fill::zeros);
  for (arma::uword r = 0; r < crosses.n_rows; ++r) {
    arma::uword a = crosses(r, 0);
    arma::uword b = crosses(r, 1);
    x(a) += val;
    x(b) += val;
  }
  return x;
}
static inline double calc_sim(const arma::vec& x, const arma::mat& G) {
  arma::uvec idx = arma::find(arma::abs(x) > 0);
  const arma::uword k = idx.n_elem;
  double q = 0.0;

  for (arma::uword a = 0; a < k; ++a) {
    const arma::uword i = idx[a];
    const double xi = x[i];
    q += xi * G(i, i) * xi;  // diagonal
    for (arma::uword b = a + 1; b < k; ++b) {
      const arma::uword j = idx[b];
      q += 2.0 * xi * G(i, j) * x[j];  // off-diagonal
    }
  }
  return q;
}

// Shared generalized constraint helpers are placed after the sampling helpers.

// --- crowding distance on a front (indices into candidate list) ---
static inline void crowding_distance(
    const std::vector<int>& front,
    const arma::vec& u,
    const arma::vec& sim,
    const std::vector<int>& cand,
    std::vector<double>& cd_out)
{
  const int m = (int)front.size();
  if (m == 0) return;

  for (int fi : front) cd_out[fi] = 0.0;
  if (m <= 2) {
    for (int fi : front) cd_out[fi] = std::numeric_limits<double>::infinity();
    return;
  }

  std::vector<int> by_u = front;
  std::sort(by_u.begin(), by_u.end(), [&](int a, int b){
    return u(cand[a]) > u(cand[b]);
  });
  std::vector<int> by_s = front;
  std::sort(by_s.begin(), by_s.end(), [&](int a, int b){
    return sim(cand[a]) < sim(cand[b]);
  });

  double umin = u(cand[by_u.back()]), umax = u(cand[by_u.front()]);
  double smin = sim(cand[by_s.front()]), smax = sim(cand[by_s.back()]);
  double urange = std::max(umax - umin, 1e-12);
  double srange = std::max(smax - smin, 1e-12);

  cd_out[by_u.front()] = cd_out[by_u.back()] = std::numeric_limits<double>::infinity();
  cd_out[by_s.front()] = cd_out[by_s.back()] = std::numeric_limits<double>::infinity();

  for (int k = 1; k < m - 1; ++k) {
    int i = by_u[k];
    if (!std::isinf(cd_out[i])) {
      double up = u(cand[by_u[k-1]]);
      double un = u(cand[by_u[k+1]]);
      cd_out[i] += (up - un) / urange;
    }
  }
  for (int k = 1; k < m - 1; ++k) {
    int i = by_s[k];
    if (!std::isinf(cd_out[i])) {
      double sp = sim(cand[by_s[k+1]]);
      double sn = sim(cand[by_s[k-1]]);
      cd_out[i] += (sp - sn) / srange;
    }
  }
}

static inline bool dominates_u_sim(double u1, double s1, double u2, double s2) {
  return (u1 >= u2 && s1 <= s2) && (u1 > u2 || s1 < s2);
}

// NOTE: this nsga2 select assumes cand is a list of indices into u/sim vectors.
// It selects nondominated-by-NSGA2 sorting + crowding. We will use it ONLY on
// the already Pareto-filtered set (so "at most by pareto" is satisfied).
static inline std::vector<int> nsga2_select_maxu_minsim(
    const arma::vec& u,
    const arma::vec& sim,
    const std::vector<int>& cand,
    int nSel)
{
  const int K = (int)cand.size();
  nSel = std::min(nSel, K);

  std::vector<int> domCount(K, 0);
  std::vector<std::vector<int>> domList(K);

  for (int p = 0; p < K; ++p) {
    const int ip = cand[p];
    for (int q = 0; q < K; ++q) if (q != p) {
      const int iq = cand[q];
      if (dominates_u_sim(u(ip), sim(ip), u(iq), sim(iq))) {
        domList[p].push_back(q);
      } else if (dominates_u_sim(u(iq), sim(iq), u(ip), sim(ip))) {
        domCount[p] += 1;
      }
    }
  }

  std::vector<std::vector<int>> fronts;
  fronts.reserve(64);

  std::vector<int> F;
  F.reserve(K);
  for (int i = 0; i < K; ++i) if (domCount[i] == 0) F.push_back(i);
  fronts.push_back(F);

  int fidx = 0;
  while (!fronts[fidx].empty()) {
    std::vector<int> Q;
    Q.reserve(K);
    for (int p : fronts[fidx]) {
      for (int q : domList[p]) {
        domCount[q] -= 1;
        if (domCount[q] == 0) Q.push_back(q);
      }
    }
    ++fidx;
    fronts.push_back(Q);
  }

  std::vector<int> selected;
  selected.reserve(nSel);
  std::vector<double> cd(K, 0.0);

  for (const auto& front : fronts) {
    if (front.empty()) break;

    crowding_distance(front, u, sim, cand, cd);

    std::vector<int> temp = front;
    std::sort(temp.begin(), temp.end(), [&](int a, int b){
      return cd[a] > cd[b];
    });

    int can_take = std::min<int>((int)temp.size(), nSel - (int)selected.size());
    for (int k = 0; k < can_take; ++k) selected.push_back(cand[temp[k]]);
    if ((int)selected.size() >= nSel) break;
  }

  return selected;
}

// map linear index -> row/col in upper triangle of n x n (row <= col)
static inline arma::uword mapRow(const arma::uword& k, const arma::uword& n) {
  return n - 2 - static_cast<arma::uword>(
      std::floor(
        std::sqrt(-8.0 * static_cast<double>(k)
                    + 4.0 * static_cast<double>(n) * (static_cast<double>(n) - 1.0) - 7.0) / 2.0
  - 0.5
      )
  );
}
static inline arma::uword mapCol(const arma::uword& row, const arma::uword& k, const arma::uword& n) {
  return k + row + 1 - n * (n - 1) / 2 + (n - row) * ((n - row) - 1) / 2;
}

// sampling without replacement: returns sorted unique indices in [0, N-1]
static inline arma::uvec sampleInt_std(
    arma::uword n,
    arma::uword N,
    std::mt19937& gen)
{
  if (n > N) n = N;

  std::vector<arma::uword> pool(N);
  std::iota(pool.begin(), pool.end(), static_cast<arma::uword>(0));

  for (arma::uword i = 0; i < n; ++i) {
    std::uniform_int_distribution<arma::uword> dis(i, N - 1);

    const arma::uword j = dis(gen);
    std::swap(pool[i], pool[j]);
  }

  std::sort(
    pool.begin(),
    pool.begin() + static_cast<std::ptrdiff_t>(n)
  );

  arma::uvec out(n);

  for (arma::uword i = 0; i < n; ++i) {
    out(i) = pool[static_cast<std::size_t>(i)];
  }

  return out;
}
// half-diallel sample (n pairs from all combinations, without replacement)
static inline arma::umat sampHalfDialComb_std(arma::uword nLevel, arma::uword n, std::mt19937& gen) {
  const arma::uword N = nLevel * (nLevel - 1) / 2;
  arma::uvec samples = sampleInt_std(std::min(n, N), N, gen);
  arma::umat out(2, samples.n_elem);
  for (arma::uword i = 0; i < samples.n_elem; ++i) {
    arma::uword r = mapRow(samples(i), nLevel);
    arma::uword c = mapCol(r, samples(i), nLevel);
    out(0, i) = r;
    out(1, i) = c;
  }
  return out;
}

// -----------------------------------------------------------------------------
// Generalized parental and cross-level constraint machinery
// Same implementation as the angle optimizer.
// -----------------------------------------------------------------------------

// Violation of a single lower/upper bound.
static inline double bound_violation(
    double value,
    double lower,
    double upper,
    double tol = 1e-12
) {
  double out = 0.0;
  
  if (value < lower - tol) {
    out += lower - value;
  }
  
  if (value > upper + tol) {
    out += value - upper;
  }
  
  return out;
}


// -----------------------------------------------------------------------------
// Individual feasibility checks
// These are necessary, but not sufficient, conditions for joint feasibility.
// -----------------------------------------------------------------------------

static inline void check_constraint_feasibility(
    const arma::umat& crosses,
    const arma::umat& fixedCrosses,
    arma::uword nCross,
    arma::uword nInd,
    const arma::vec& minContr,
    const arma::vec& maxContr,
    const arma::mat& A,
    const arma::mat& Afixed,
    const arma::vec& lower,
    const arma::vec& upper,
    double tol = 1e-12
) {
  
  const arma::uword nFixed = fixedCrosses.n_rows;
  const arma::uword nVar = nCross - nFixed;
  const arma::uword potCross = crosses.n_rows;
  const arma::uword K = A.n_cols;
  
  // Check dimensions.
  if (minContr.n_elem != nInd || maxContr.n_elem != nInd) {
    Rcpp::stop(
      "minContr and maxContr must have one entry per parent."
    );
  }
  
  if (
      A.n_rows != potCross ||
        Afixed.n_rows != nFixed ||
        Afixed.n_cols != K ||
        lower.n_elem != K ||
        upper.n_elem != K
  ) {
    Rcpp::stop(
      "Additional constraint dimensions are inconsistent."
    );
  }
  
  // Validate parental contribution limits.
  double sumLower = 0.0;
  double sumUpper = 0.0;
  
  for (arma::uword p = 0; p < nInd; ++p) {
    
    if (
        !std::isfinite(minContr(p)) ||
          minContr(p) < 0 ||
          minContr(p) > 1 ||
          std::isnan(maxContr(p)) ||
          maxContr(p) < 0 ||
          minContr(p) > maxContr(p) + tol
    ) {
      Rcpp::stop(
        "Invalid parental contribution bounds for parent %u.",
        static_cast<unsigned>(p + 1)
      );
    }
    
    sumLower += minContr(p);
    sumUpper += std::min(1.0, maxContr(p));
  }
  
  if (sumLower > 1.0 + tol || sumUpper < 1.0 - tol) {
    Rcpp::stop(
      "Parental contribution bounds cannot sum to one."
    );
  }
  
  // Validate additional constraint limits.
  for (arma::uword j = 0; j < K; ++j) {
    
    if (
        !std::isfinite(lower(j)) ||
          lower(j) < 0 ||
          lower(j) > 1 ||
          std::isnan(upper(j)) ||
          upper(j) < 0 ||
          lower(j) > upper(j) + tol
    ) {
      Rcpp::stop(
        "Invalid bounds for additional constraint %u.",
        static_cast<unsigned>(j + 1)
      );
    }
  }
  
  // Validate the additional cross scores.
  // All values must be frequencies or proportions between 0 and 1.
  for (arma::uword i = 0; i < potCross; ++i) {
    for (arma::uword j = 0; j < K; ++j) {
      
      if (
          !std::isfinite(A(i, j)) ||
            A(i, j) < 0 ||
            A(i, j) > 1
      ) {
        Rcpp::stop(
          "A must contain finite values in [0,1]."
        );
      }
    }
  }
  
  for (arma::uword i = 0; i < nFixed; ++i) {
    for (arma::uword j = 0; j < K; ++j) {
      
      if (
          !std::isfinite(Afixed(i, j)) ||
            Afixed(i, j) < 0 ||
            Afixed(i, j) > 1
      ) {
        Rcpp::stop(
          "Afixed must contain finite values in [0,1]."
        );
      }
    }
  }
  

  // For each candidate cross, a parent can appear 0, 1, or 2 times.
  // Count these occurrences without allocating a dense parent matrix.
  
  std::vector<arma::uword> countOne(nInd, 0);
  std::vector<arma::uword> countTwo(nInd, 0);
  
  for (arma::uword i = 0; i < potCross; ++i) {
    
    const arma::uword a = crosses(i, 0);
    const arma::uword b = crosses(i, 1);
    
    if (a == b) {
      ++countTwo[a];
    } else {
      ++countOne[a];
      ++countOne[b];
    }
  }
  
  // Fixed parental contributions.
  arma::vec fixedParent =
    calcContr(fixedCrosses, nInd, nCross);
  
  for (arma::uword p = 0; p < nInd; ++p) {
    
    const arma::uword n0 =
      potCross - countOne[p] - countTwo[p];
    
    // Minimum possible contribution:
    // select crosses with 0 appearances first,
    // followed by 1, followed by 2.
    
    const arma::uword take0 =
      std::min(nVar, n0);
    
    const arma::uword take1 =
      std::min(nVar - take0, countOne[p]);
    
    const arma::uword take2 =
      nVar - take0 - take1;
    
    const double minPossible =
      fixedParent(p) +
      (
          static_cast<double>(take1) +
            2.0 * take2
      ) / (2.0 * nCross);
    
    // Maximum possible contribution:
    // select crosses with 2 appearances first,
    // followed by 1, followed by 0.
    
    const arma::uword take2max =
      std::min(nVar, countTwo[p]);
    
    const arma::uword take1max =
      std::min(nVar - take2max, countOne[p]);
    
    const double maxPossible =
      fixedParent(p) +
      (
          2.0 * take2max +
            static_cast<double>(take1max)
      ) / (2.0 * nCross);
    
    if (
        minContr(p) > maxPossible + tol ||
          maxContr(p) < minPossible - tol
    ) {
      Rcpp::stop(
        "Parent %u impossible individually: "
        "achievable [%.6f, %.6f], requested [%.6f, %.6f].",
        static_cast<unsigned>(p + 1),
        minPossible,
        maxPossible,
        minContr(p),
        maxContr(p)
      );
    }
  }
  
  // ---------------------------------------------------------------------------
  // Feasibility of additional constraints
  // ---------------------------------------------------------------------------
  
  for (arma::uword j = 0; j < K; ++j) {
    
    // Contribution from fixed crosses.
    double fixedSum = 0.0;
    
    for (arma::uword f = 0; f < nFixed; ++f) {
      fixedSum += Afixed(f, j);
    }
    
    // Collect scores for the current constraint.
    std::vector<double> scores(
        static_cast<size_t>(potCross)
    );
    
    for (arma::uword i = 0; i < potCross; ++i) {
      scores[i] = A(i, j);
    }
    
    double minSum = 0.0;
    double maxSum = 0.0;
    
    if (nVar == potCross) {
      
      // All candidate crosses must be included.
      minSum = std::accumulate(
        scores.begin(),
        scores.end(),
        0.0
      );
      
      maxSum = minSum;
      
    } else if (nVar > 0) {
      
      // Select the nVar smallest scores.
      std::nth_element(
        scores.begin(),
        scores.begin() + nVar,
        scores.end()
      );
      
      minSum = std::accumulate(
        scores.begin(),
        scores.begin() + nVar,
        0.0
      );
      
      // Select the nVar largest scores.
      std::nth_element(
        scores.begin(),
        scores.end() - nVar,
        scores.end()
      );
      
      maxSum = std::accumulate(
        scores.end() - nVar,
        scores.end(),
        0.0
      );
    }
    
    const double minPossible =
      (fixedSum + minSum) / nCross;
    
    const double maxPossible =
      (fixedSum + maxSum) / nCross;
    
    if (
        lower(j) > maxPossible + tol ||
          upper(j) < minPossible - tol
    ) {
      Rcpp::stop(
        "Additional constraint %u impossible individually: "
        "achievable [%.6f, %.6f], requested [%.6f, %.6f].",
        static_cast<unsigned>(j + 1),
        minPossible,
        maxPossible,
        lower(j),
        upper(j)
      );
    }
  }
}


struct ConstraintState {
  
  arma::vec parents;  // Parental contributions
  arma::vec extra;    // QTL frequencies and cross-type proportions
  
  double violation = 0.0;
};



struct ConstraintSystem {
  
  const arma::umat& crosses;
  const arma::mat& scores;      // K x potCross
  
  const arma::vec& minContr;
  const arma::vec& maxContr;
  
  const arma::vec& lower;
  const arma::vec& upper;
  
  const arma::vec& fixedParents;
  const arma::vec& fixedExtra;
  
  arma::uword nCross;
  arma::uword nInd;
  
  double tol;
  
  
  // Check whether any constraint is active.
  bool active() const {
    
    for (arma::uword p = 0; p < nInd; ++p) {
      
      if (
          minContr(p) > tol ||
            maxContr(p) < 1.0 - tol
      ) {
        return true;
      }
    }
    
    for (arma::uword j = 0; j < lower.n_elem; ++j) {
      
      if (
          lower(j) > tol ||
            upper(j) < 1.0 - tol
      ) {
        return true;
      }
    }
    
    return false;
  }
  
  
  // Calculate total violation across all constraints.
  double full_violation(
      const ConstraintState& s
  ) const {
    
    double v = 0.0;
    
    for (arma::uword p = 0; p < nInd; ++p) {
      
      v += bound_violation(
        s.parents(p),
        minContr(p),
        maxContr(p),
        tol
      );
    }
    
    for (arma::uword j = 0; j < lower.n_elem; ++j) {
      
      v += bound_violation(
        s.extra(j),
        lower(j),
        upper(j),
        tol
      );
    }
    
    return v;
  }
  
  
  // Initialize constraints for a complete variable-cross plan.
  // Fixed crosses are already included through fixedParents/fixedExtra.
  ConstraintState initialize(
      const arma::uvec& plan
  ) const {
    
    ConstraintState s;
    
    s.parents = fixedParents;
    s.extra = fixedExtra;
    
    const double parentStep =
      1.0 / (2.0 * nCross);
    
    const double extraStep =
      1.0 / nCross;
    
    for (arma::uword k = 0; k < plan.n_elem; ++k) {
      
      const arma::uword c = plan(k);
      
      s.parents(crosses(c, 0)) += parentStep;
      s.parents(crosses(c, 1)) += parentStep;
      
      for (arma::uword j = 0; j < lower.n_elem; ++j) {
        s.extra(j) += scores(j, c) * extraStep;
      }
    }
    
    s.violation = full_violation(s);
    
    return s;
  }
  
  
  // Evaluate a proposed single-cross swap without changing the state.
  double after_swap(
      const ConstraintState& s,
      arma::uword oldCross,
      arma::uword newCross
  ) const {
    
    const arma::uword oa = crosses(oldCross, 0);
    const arma::uword ob = crosses(oldCross, 1);
    const arma::uword na = crosses(newCross, 0);
    const arma::uword nb = crosses(newCross, 1);
    
    arma::uword affected[4] = {oa, ob, na, nb};
    
    const double step =
      1.0 / (2.0 * nCross);
    
    double result = s.violation;
    
    // Parent contributions: inspect at most four distinct parents.
    for (int i = 0; i < 4; ++i) {
      
      bool seen = false;
      
      for (int k = 0; k < i; ++k) {
        if (affected[k] == affected[i]) {
          seen = true;
          break;
        }
      }
      
      if (seen) continue;
      
      const arma::uword p = affected[i];
      
      const double delta = step * (
        static_cast<int>(na == p) +
          static_cast<int>(nb == p) -
          static_cast<int>(oa == p) -
          static_cast<int>(ob == p)
      );
      
      result +=
        bound_violation(
          s.parents(p) + delta,
          minContr(p),
          maxContr(p),
          tol
        ) -
          bound_violation(
            s.parents(p),
            minContr(p),
            maxContr(p),
            tol
          );
    }
    
    // Additional constraints: evaluate changed frequencies.
    for (arma::uword j = 0; j < lower.n_elem; ++j) {
      
      const double next =
        s.extra(j) +
        (scores(j, newCross) - scores(j, oldCross)) / nCross;
      
      result +=
        bound_violation(
          next,
          lower(j),
          upper(j),
          tol
        ) -
          bound_violation(
            s.extra(j),
            lower(j),
            upper(j),
            tol
          );
    }
    
    return std::max(0.0, result);
  }
  
  
  // Evaluate two simultaneous replacements.
  double after_double_swap(
      const ConstraintState& s,
      arma::uword old1,
      arma::uword new1,
      arma::uword old2,
      arma::uword new2
  ) const {
    
    const arma::uword affected[8] = {
      crosses(old1, 0),
      crosses(old1, 1),
      crosses(new1, 0),
      crosses(new1, 1),
      crosses(old2, 0),
      crosses(old2, 1),
      crosses(new2, 0),
      crosses(new2, 1)
    };
    
    const double step =
      1.0 / (2.0 * nCross);
    
    double result = s.violation;
    
    for (int i = 0; i < 8; ++i) {
      
      bool seen = false;
      
      for (int k = 0; k < i; ++k) {
        if (affected[k] == affected[i]) {
          seen = true;
          break;
        }
      }
      
      if (seen) continue;
      
      const arma::uword p = affected[i];
      
      const double delta = step * (
        static_cast<int>(affected[2] == p) +
          static_cast<int>(affected[3] == p) +
          static_cast<int>(affected[6] == p) +
          static_cast<int>(affected[7] == p) -
          static_cast<int>(affected[0] == p) -
          static_cast<int>(affected[1] == p) -
          static_cast<int>(affected[4] == p) -
          static_cast<int>(affected[5] == p)
      );
      
      result +=
        bound_violation(
          s.parents(p) + delta,
          minContr(p),
          maxContr(p),
          tol
        ) -
          bound_violation(
            s.parents(p),
            minContr(p),
            maxContr(p),
            tol
          );
    }
    
    for (arma::uword j = 0; j < lower.n_elem; ++j) {
      
      const double next =
        s.extra(j) +
        (
            scores(j, new1) - scores(j, old1) +
              scores(j, new2) - scores(j, old2)
        ) / nCross;
      
      result +=
        bound_violation(
          next,
          lower(j),
          upper(j),
          tol
        ) -
          bound_violation(
            s.extra(j),
            lower(j),
            upper(j),
            tol
          );
    }
    
    return std::max(0.0, result);
  }
  
  
  // Apply an accepted single-cross replacement.
  void apply_swap(
      ConstraintState& s,
      arma::uword oldCross,
      arma::uword newCross,
      double newViolation
  ) const {
    
    const double step =
      1.0 / (2.0 * nCross);
    
    s.parents(crosses(oldCross, 0)) -= step;
    s.parents(crosses(oldCross, 1)) -= step;
    
    s.parents(crosses(newCross, 0)) += step;
    s.parents(crosses(newCross, 1)) += step;
    
    for (arma::uword j = 0; j < lower.n_elem; ++j) {
      
      s.extra(j) +=
        (scores(j, newCross) - scores(j, oldCross)) / nCross;
    }
    
    s.violation = newViolation;
  }
};

// -----------------------------------------------------------------------------
// Repair a single mating plan
// -----------------------------------------------------------------------------

static inline bool repair_one_plan(
    arma::uvec& plan,
    const ConstraintSystem& cs,
    std::mt19937& gen,
    int maxIter = 2000,
    int maxCandidateTries = 3000
) {
  
  ConstraintState state = cs.initialize(plan);
  
  if (state.violation <= cs.tol) {
    return true;
  }
  
  if (
      plan.n_elem == 0 ||
        cs.crosses.n_rows <= plan.n_elem
  ) {
    return false;
  }
  
  // Track selected crosses to enforce uniqueness.
  std::unordered_set<arma::uword> selected(
      plan.begin(),
      plan.end()
  );
  
  std::uniform_int_distribution<arma::uword> posDist(
      0,
      plan.n_elem - 1
  );
  
  std::uniform_int_distribution<arma::uword> crossDist(
      0,
      cs.crosses.n_rows - 1
  );
  
  for (
      int iter = 0;
      iter < maxIter && state.violation > cs.tol;
      ++iter
  ) {
    
    bool improved = false;
    
    arma::uword chosenPos = 0;
    arma::uword chosenNew = 0;
    
    double chosenV = state.violation;
    
    // -------------------------------------------------------------------------
    // First attempt: random single-cross replacements.
    // -------------------------------------------------------------------------
    
    for (int t = 0; t < maxCandidateTries; ++t) {
      
      const arma::uword pos = posDist(gen);
      const arma::uword newCross = crossDist(gen);
      
      // Do not introduce duplicate crosses.
      if (selected.count(newCross)) {
        continue;
      }
      
      const double newV =
        cs.after_swap(
          state,
          plan(pos),
          newCross
        );
      
      if (newV < state.violation - cs.tol) {
        
        chosenPos = pos;
        chosenNew = newCross;
        chosenV = newV;
        
        improved = true;
        break;
      }
    }
    
    // -------------------------------------------------------------------------
    // Second attempt: exhaustive single-cross replacements.
    // Only used for manageable numbers of combinations.
    // -------------------------------------------------------------------------
    
    if (
        !improved &&
          static_cast<double>(plan.n_elem) *
          cs.crosses.n_rows <= 200000.0
    ) {
      
      for (
          arma::uword c = 0;
          c < cs.crosses.n_rows && !improved;
          ++c
      ) {
        
        if (selected.count(c)) {
          continue;
        }
        
        for (arma::uword pos = 0; pos < plan.n_elem; ++pos) {
          
          const double newV =
            cs.after_swap(
              state,
              plan(pos),
              c
            );
          
          if (newV < state.violation - cs.tol) {
            
            chosenPos = pos;
            chosenNew = c;
            chosenV = newV;
            
            improved = true;
            break;
          }
        }
      }
    }
    
    // -------------------------------------------------------------------------
    // Third attempt: simultaneous two-cross replacement.
    // This can overcome local minima of single-cross repair.
    // -------------------------------------------------------------------------
    
    if (
        !improved &&
          plan.n_elem >= 2 &&
          cs.crosses.n_rows >= plan.n_elem + 2
    ) {
      
      for (int trial = 0; trial < 2500; ++trial) {
        
        const arma::uword pos1 = posDist(gen);
        const arma::uword pos2 = posDist(gen);
        
        if (pos1 == pos2) {
          continue;
        }
        
        const arma::uword new1 = crossDist(gen);
        const arma::uword new2 = crossDist(gen);
        
        if (
            new1 == new2 ||
              selected.count(new1) ||
              selected.count(new2)
        ) {
          continue;
        }
        
        const double newV =
          cs.after_double_swap(
            state,
            plan(pos1),
            new1,
            plan(pos2),
            new2
          );
        
        if (newV < state.violation - cs.tol) {
          
          const arma::uword old1 = plan(pos1);
          const arma::uword old2 = plan(pos2);
          
          // Apply the first replacement.
          const double midV =
            cs.after_swap(
              state,
              old1,
              new1
            );
          
          cs.apply_swap(
            state,
            old1,
            new1,
            midV
          );
          
          // Apply the second replacement.
          cs.apply_swap(
            state,
            old2,
            new2,
            newV
          );
          
          selected.erase(old1);
          selected.erase(old2);
          
          selected.insert(new1);
          selected.insert(new2);
          
          plan(pos1) = new1;
          plan(pos2) = new2;
          
          improved = true;
          break;
        }
      }
      
      if (improved) {
        
        if ((iter + 1) % 64 == 0) {
          state.violation = cs.full_violation(state);
        }
        
        continue;
      }
    }
    
    // No improving replacement was found.
    if (!improved) {
      break;
    }
    
    // -------------------------------------------------------------------------
    // Apply the accepted single-cross replacement.
    // -------------------------------------------------------------------------
    
    const arma::uword oldCross = plan(chosenPos);
    
    cs.apply_swap(
      state,
      oldCross,
      chosenNew,
      chosenV
    );
    
    selected.erase(oldCross);
    selected.insert(chosenNew);
    
    plan(chosenPos) = chosenNew;
    
    // Periodically recalculate violations to limit numerical drift.
    if ((iter + 1) % 64 == 0) {
      state.violation = cs.full_violation(state);
    }
  }
  
  // Recalculate the entire constraint state to verify feasibility.
  return cs.initialize(plan).violation <= cs.tol;
}

// -----------------------------------------------------------------------------
// Diversify a feasible mating plan
// -----------------------------------------------------------------------------

static inline void diversify_feasible_plan(
    arma::uvec& plan,
    const ConstraintSystem& cs,
    std::mt19937& gen,
    int maxAccepted = 3
) {
  
  if (
      plan.n_elem == 0 ||
        cs.crosses.n_rows <= plan.n_elem
  ) {
    return;
  }
  
  ConstraintState s = cs.initialize(plan);
  
  std::unordered_set<arma::uword> selected(
      plan.begin(),
      plan.end()
  );
  
  std::uniform_int_distribution<arma::uword> posDist(
      0,
      plan.n_elem - 1
  );
  
  std::uniform_int_distribution<arma::uword> crossDist(
      0,
      cs.crosses.n_rows - 1
  );
  
  int accepted = 0;
  
  // Random feasible single-cross replacements.
  for (
      int t = 0;
      t < 200 && accepted < maxAccepted;
      ++t
  ) {
    
    const arma::uword pos = posDist(gen);
    const arma::uword next = crossDist(gen);
    
    if (selected.count(next)) {
      continue;
    }
    
    const arma::uword old = plan(pos);
    
    const double v =
      cs.after_swap(
        s,
        old,
        next
      );
    
    if (v <= cs.tol) {
      
      cs.apply_swap(
        s,
        old,
        next,
        v
      );
      
      selected.erase(old);
      selected.insert(next);
      
      plan(pos) = next;
      
      ++accepted;
    }
  }
  
  // If no feasible single swap was possible,
  // attempt simultaneous two-cross replacements.
  
  if (
      accepted == 0 &&
        plan.n_elem >= 2 &&
        cs.crosses.n_rows >= plan.n_elem + 2
  ) {
    
    for (int t = 0; t < 500; ++t) {
      
      const arma::uword pos1 = posDist(gen);
      const arma::uword pos2 = posDist(gen);
      
      if (pos1 == pos2) {
        continue;
      }
      
      const arma::uword new1 = crossDist(gen);
      const arma::uword new2 = crossDist(gen);
      
      if (
          new1 == new2 ||
            selected.count(new1) ||
            selected.count(new2)
      ) {
        continue;
      }
      
      const arma::uword old1 = plan(pos1);
      const arma::uword old2 = plan(pos2);
      
      const double newV =
        cs.after_double_swap(
          s,
          old1,
          new1,
          old2,
          new2
        );
      
      if (newV <= cs.tol) {
        
        const double midV =
          cs.after_swap(
            s,
            old1,
            new1
          );
        
        cs.apply_swap(
          s,
          old1,
          new1,
          midV
        );
        
        cs.apply_swap(
          s,
          old2,
          new2,
          newV
        );
        
        selected.erase(old1);
        selected.erase(old2);
        
        selected.insert(new1);
        selected.insert(new2);
        
        plan(pos1) = new1;
        plan(pos2) = new2;
        
        break;
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Construct a feasible high-gain reference plan
// -----------------------------------------------------------------------------

static inline arma::uvec build_constrained_umax_plan(
    const arma::vec& u,
    const arma::uvec& uOrder,
    const ConstraintSystem& cs,
    arma::uword nVar,
    uint64_t base_seed
) {
  
  // Start with the nVar crosses with highest gain.
  arma::uvec plan = uOrder.head(nVar);
  
  // This is already the unconstrained gain maximum.
  if (
      !cs.active() ||
        cs.initialize(plan).violation <= cs.tol
  ) {
    return plan;
  }
  
  std::mt19937 gen(
      static_cast<uint32_t>(
        base_seed ^ 0xC44D61A7ULL
      )
  );
  
  // Attempt to repair the highest-gain plan.
  bool found =
    repair_one_plan(
      plan,
      cs,
      gen
    );
  
  // Try alternative starting plans if repair failed.
  for (
      int attempt = 0;
      !found && attempt < 20;
      ++attempt
  ) {
    
    plan = sampleInt_std(
      nVar,
      cs.crosses.n_rows,
      gen
    );
    
    found =
      repair_one_plan(
        plan,
        cs,
        gen
      );
  }
  
  if (!found) {
    
    Rcpp::stop(
      "Cannot construct a jointly feasible mating plan. "
      "Each bound passed the individual feasibility checks, "
      "but the constraints may conflict or the swap heuristic "
      "may be stuck."
    );
  }
  
  // ---------------------------------------------------------------------------
  // Improve gain while keeping all constraints feasible.
  // ---------------------------------------------------------------------------
  
  std::unordered_set<arma::uword> selected(
      plan.begin(),
      plan.end()
  );
  
  ConstraintState state = cs.initialize(plan);
  
  for (
      arma::uword iteration = 0;
      iteration < 10 * nVar;
      ++iteration
  ) {
    
    bool improved = false;
    
    // Prefer removing selected crosses with low gain.
    std::vector<arma::uword> removeOrder(plan.n_elem);
    
    std::iota(
      removeOrder.begin(),
      removeOrder.end(),
      arma::uword(0)
    );
    
    std::sort(
      removeOrder.begin(),
      removeOrder.end(),
      [&](arma::uword a, arma::uword b) {
        return u(plan(a)) < u(plan(b));
      }
    );
    
    // Candidate insertions are ordered by decreasing gain.
    for (
        arma::uword k = 0;
        k < uOrder.n_elem && !improved;
        ++k
    ) {
      
      const arma::uword next = uOrder(k);
      
      if (selected.count(next)) {
        continue;
      }
      
      for (arma::uword pos : removeOrder) {
        
        const arma::uword old = plan(pos);
        
        if (u(next) <= u(old) + cs.tol) {
          break;
        }
        
        const double v =
          cs.after_swap(
            state,
            old,
            next
          );
        
        if (v <= cs.tol) {
          
          cs.apply_swap(
            state,
            old,
            next,
            v
          );
          
          selected.erase(old);
          selected.insert(next);
          
          plan(pos) = next;
          
          improved = true;
          break;
        }
      }
    }
    
    if (!improved) {
      break;
    }
  }
  
  if (cs.initialize(plan).violation > cs.tol) {
    Rcpp::stop(
      "Internal error: maximum-gain plan violates constraints."
    );
  }
  
  return plan;
}

// -----------------------------------------------------------------------------
// Repair all mating plans in a population
// -----------------------------------------------------------------------------

static inline void repair_constraint_plans(
    arma::umat& Plans,
    const ConstraintSystem& cs,
    const arma::uvec& feasibleSeed,
    uint64_t base_seed
) {
  
  if (!cs.active() || Plans.n_cols == 0) {
    return;
  }
  
  std::vector<int> usedFallback(
      static_cast<size_t>(Plans.n_cols),
      0
  );
  
  ct_parallel_for(
    0,
    static_cast<int>(Plans.n_cols),
    [&](int ii) {
      
      const arma::uword planIndex =
        static_cast<arma::uword>(ii);
      
      std::mt19937 gen(
          static_cast<uint32_t>(
            base_seed ^
              (
                  0xA13B77C19E3779B9ULL +
                    static_cast<uint64_t>(ii)
              )
          )
      );
      
      arma::uvec plan =
        Plans.col(planIndex);
      
      // Repair the current plan.
      if (!repair_one_plan(plan, cs, gen)) {
        
        // If repair fails, start from a known feasible plan.
        plan = feasibleSeed;
        
        // Attempt to introduce feasible variation.
        diversify_feasible_plan(
          plan,
          cs,
          gen
        );
        
        usedFallback[static_cast<size_t>(ii)] = 1;
      }
      
      Plans.col(planIndex) = plan;
    }
  );
  
  // Report how many individuals required feasible-seed fallback.
  const size_t nFallback =
    std::accumulate(
      usedFallback.begin(),
      usedFallback.end(),
      size_t(0)
    );
  
  if (nFallback > 0) {
    
    // Rcpp::Rcout
    // << "Constraint repair: "
    // << nFallback
    // << " individual(s) replaced by feasible-seed variants."
    // << std::endl;
  }
}



// uniform set-based crossover
static inline arma::uvec mate_uniform(
    const arma::uvec& parent1,
    const arma::uvec& parent2,
    std::mt19937& gen,
    arma::uword potCross)
{
  const arma::uword n = parent1.n_elem;

  if (parent2.n_elem != n) {
    Rcpp::stop(
      "Internal error: mating parents have different lengths."
    );
  }

  if (n == 0) {
    return arma::uvec();
  }

  /*
   * Shared crosses are inherited automatically.
   * differing contains crosses found in only one parent.
   */
  std::vector<arma::uword> child;
  std::vector<arma::uword> differing;

  child.reserve(static_cast<size_t>(n));
  differing.reserve(static_cast<size_t>(2 * n));

  /*
   * Process parent 1:
   *
   * - if a cross is also in parent 2, add it to the child;
   * - otherwise add it to the pool of differing parental crosses.
   */
  for (arma::uword i = 0; i < n; ++i) {
    const arma::uword gene = parent1(i);

    bool shared = false;

    for (arma::uword j = 0; j < n; ++j) {
      if (parent2(j) == gene) {
        shared = true;
        break;
      }
    }

    if (shared) {
      child.push_back(gene);
    } else {
      differing.push_back(gene);
    }
  }

  /*
   * Add crosses unique to parent 2.
   */
  for (arma::uword i = 0; i < n; ++i) {
    const arma::uword gene = parent2(i);

    bool shared = false;

    for (arma::uword j = 0; j < n; ++j) {
      if (parent1(j) == gene) {
        shared = true;
        break;
      }
    }

    if (!shared) {
      differing.push_back(gene);
    }
  }

  /*
   * Randomize the differing parental crosses and take enough
   * to fill the offspring.
   */
  std::shuffle(
    differing.begin(),
    differing.end(),
    gen
  );

  const arma::uword needed =
    n - static_cast<arma::uword>(child.size());

  if (differing.size() < static_cast<size_t>(needed)) {
    Rcpp::stop(
      "Internal error: parental union is too small to construct offspring."
    );
  }

  for (arma::uword i = 0; i < needed; ++i) {
    child.push_back(
      differing[static_cast<size_t>(i)]
    );
  }

  /*
   * The plan is conceptually unordered, but shuffling prevents
   * shared crosses from always occupying the first positions.
   */
  std::shuffle(
    child.begin(),
    child.end(),
    gen
  );

  arma::uvec offspring(n);

  for (arma::uword i = 0; i < n; ++i) {
    offspring(i) =
      child[static_cast<size_t>(i)];
  }

  return offspring;
}

// mutation: replace nMutate positions with genes not in plan
static inline arma::uvec mutate(
    const arma::uvec& crosses,
    arma::uword nMutate,
    arma::uword potCross,
    std::mt19937& gen)
{
  const arma::uword n = crosses.n_elem;

  if (n == 0 || potCross <= n) {
    return crosses;
  }

  const arma::uword k =
    std::min<arma::uword>(nMutate, n);

  if (k == 0) {
    return crosses;
  }

  arma::uvec out = crosses;

  /*
   * This samples positions from the plan, so N is only nVar,
   * not potCross.
   */
  arma::uvec positions =
    sampleInt_std(k, n, gen);

  std::uniform_int_distribution<arma::uword> gene_dist(
      0,
      potCross - 1
  );

  /*
   * Check only the current plan for duplicates.
   */
  const auto used_elsewhere =
    [&](arma::uword gene, arma::uword skip) -> bool {
      for (arma::uword j = 0; j < n; ++j) {
        if (j != skip && out(j) == gene) {
          return true;
        }
      }

      return false;
    };

    for (arma::uword t = 0; t < positions.n_elem; ++t) {
      const arma::uword position = positions(t);
      const arma::uword old_gene = out(position);

      arma::uword new_gene;

      do {
        new_gene = gene_dist(gen);
      } while (
          new_gene == old_gene ||
            used_elsewhere(new_gene, position)
      );

      out(position) = new_gene;
    }

    return out;
}

// Build a stable hash key for a plan (sorted indices)
static inline std::string plan_key_sorted(const arma::uvec& plan_sorted) {
  std::ostringstream oss;
  for (arma::uword i = 0; i < plan_sorted.n_elem; ++i) {
    if (i) oss << ',';
    oss << plan_sorted(i);
  }
  return oss.str();
}

// Ensure uniqueness by resampling duplicates (cheap and effective)
// Preserve feasibility while attempting to remove duplicate mating plans.
// Some constraints may allow only very few distinct feasible plans; in that
// case duplicates are retained rather than introducing infeasible plans.
static inline void enforce_unique_population(
    arma::umat& Progeny,
    const ConstraintSystem& cs,
    const arma::uvec& feasibleSeed,
    uint64_t base_seed)
{
  std::unordered_set<std::string> seen;
  seen.reserve(static_cast<size_t>(Progeny.n_cols) * 2);

  const arma::uword nVar = Progeny.n_rows;
  const arma::uword potCross = cs.crosses.n_rows;

  // When every candidate is selected (or no variable cross exists),
  // a second distinct mating plan cannot be constructed.
  if (nVar == 0 || nVar == potCross) return;

  for (arma::uword i = 0; i < Progeny.n_cols; ++i) {
    arma::uvec plan = arma::conv_to<arma::uvec>::from(Progeny.col(i));

    if (cs.initialize(plan).violation > 1e-10) {
      Rcpp::stop("Internal error: nonfeasible plan in uniqueness check.");
    }

    std::string key = plan_key_sorted(arma::sort(plan));
    if (seen.insert(key).second) continue;

    std::mt19937 gen(static_cast<uint32_t>(
        base_seed ^ (0x8F3A2B7C1D9E6A5BULL +
                     static_cast<uint64_t>(i) * 1315423911ULL)));

    bool found = false;
    for (int attempt = 0; attempt < 48 && !found; ++attempt) {
      arma::uvec candidate = (attempt % 4 == 3) ? feasibleSeed : plan;

      if (cs.active()) {
        // Every accepted change preserves all constraint bounds.
        diversify_feasible_plan(candidate, cs, gen, 1 + (attempt % 3));
      } else if (nVar > 0 && potCross > nVar) {
        candidate = sampleInt_std(nVar, potCross, gen);
      }

      if (cs.initialize(candidate).violation > 1e-10) continue;

      const std::string candidateKey =
          plan_key_sorted(arma::sort(candidate));
      if (seen.insert(candidateKey).second) {
        Progeny.col(i) = candidate;
        found = true;
      }
    }
    // If no distinct feasible plan is found, keep the original duplicate.
  }
}

// Evaluate one plan -> u/sim
static inline void eval_plan(
    const arma::uvec& idx,              // indices into Crosses_mat rows
    const arma::vec& u,                 // per-cross utility (length potCross)
    double ufixedSum,
    arma::uword nCross,                 // total crosses (fixed + variable)
    const arma::mat& G,                 // nInd x nInd
    const arma::umat& Crosses_mat,      // potCross x 2 (0-based parents)
    arma::uword nInd,
    const arma::vec& xfixed,            // length nInd
    double& uOut,
    double& simOut)
{
  arma::vec x = calcContr(Crosses_mat.rows(idx), nInd,nCross);
  simOut = calc_sim(x + xfixed, G);
  uOut   = (arma::accu(u.elem(idx)) + ufixedSum) / static_cast<double>(nCross);
}


// Pareto-filter archive (nondominated only)  [maximize U, minimize S]
static inline void pareto_filter_archive(
    std::vector<arma::uvec>& plans,
    std::vector<double>& U,
    std::vector<double>& S)
{
  const size_t N = U.size();
  if (N == 0) return;

  // Do not keep repeated copies of the same feasible mating plan.
  std::unordered_set<std::string> seenPlans;
  seenPlans.reserve(N * 2);
  std::vector<char> keep(N, 1);
  for (size_t i = 0; i < N; ++i) {
    const std::string key = plan_key_sorted(arma::sort(plans[i]));
    if (!seenPlans.insert(key).second) keep[i] = 0;
  }
  for (size_t i = 0; i < N; ++i) {
    if (!keep[i]) continue;
    for (size_t j = 0; j < N; ++j) {
      if (i == j || !keep[j]) continue;
      if (dominates_u_sim(U[j], S[j], U[i], S[i])) { keep[i] = 0; break; }
    }
  }

  std::vector<arma::uvec> plans2;
  std::vector<double> U2, S2;
  plans2.reserve(N); U2.reserve(N); S2.reserve(N);

  for (size_t i = 0; i < N; ++i) if (keep[i]) {
    plans2.push_back(plans[i]);
    U2.push_back(U[i]);
    S2.push_back(S[i]);
  }

  plans.swap(plans2);
  U.swap(U2);
  S.swap(S2);
}


static inline void cap_archive_by_arc_length(
    std::vector<arma::uvec>& plans,
    std::vector<double>& U,
    std::vector<double>& S,
    size_t cap)
{
  const size_t N = U.size();
  if (N <= cap) return;



  // Sort archive points along the Pareto front by similarity
  std::vector<size_t> ord(N);
  for (size_t i = 0; i < N; ++i) ord[i] = i;

  std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {
    if (S[a] != S[b]) return S[a] < S[b];
    return U[a] > U[b];
  });

  // Normalize objectives before measuring distances
  double umin = *std::min_element(U.begin(), U.end());
  double umax = *std::max_element(U.begin(), U.end());
  double smin = *std::min_element(S.begin(), S.end());
  double smax = *std::max_element(S.begin(), S.end());

  double ur = std::max(umax - umin, 1e-12);
  double sr = std::max(smax - smin, 1e-12);

  // Cumulative arc length along the sorted Pareto front
  std::vector<double> arc(N, 0.0);

  for (size_t k = 1; k < N; ++k) {
    size_t i0 = ord[k - 1];
    size_t i1 = ord[k];

    double du = (U[i1] - U[i0]) / ur;
    double ds = (S[i1] - S[i0]) / sr;

    arc[k] = arc[k - 1] + std::sqrt(du * du + ds * ds);
  }

  double total_arc = arc[N - 1];

  // Fallback if all points are effectively identical
  if (total_arc <= 1e-12) {
    std::vector<arma::uvec> plans2;
    std::vector<double> U2, S2;

    plans2.reserve(cap);
    U2.reserve(cap);
    S2.reserve(cap);

    for (size_t k = 0; k < cap; ++k) {
      size_t i = ord[k];
      plans2.push_back(plans[i]);
      U2.push_back(U[i]);
      S2.push_back(S[i]);
    }

    plans.swap(plans2);
    U.swap(U2);
    S.swap(S2);
    return;
  }

  // Select points closest to equally spaced arc-length targets
  std::vector<char> keep(N, 0);
  size_t last_pos = 0;

  for (size_t m = 0; m < cap; ++m) {
    double target = total_arc * static_cast<double>(m) /
      static_cast<double>(cap - 1);

    while (last_pos + 1 < N && arc[last_pos + 1] < target) {
      ++last_pos;
    }

    size_t best_pos = last_pos;

    if (last_pos + 1 < N) {
      double d0 = std::abs(arc[last_pos] - target);
      double d1 = std::abs(arc[last_pos + 1] - target);
      if (d1 < d0) best_pos = last_pos + 1;
    }

    keep[ord[best_pos]] = 1;
  }

  // If duplicates caused fewer than cap selected points, fill gaps
  size_t nKeep = 0;
  for (char k : keep) if (k) ++nKeep;

  if (nKeep < cap) {
    std::vector<size_t> pos_order(N);
    for (size_t k = 0; k < N; ++k) pos_order[k] = k;

    std::sort(pos_order.begin(), pos_order.end(), [&](size_t a, size_t b) {
      double da = std::numeric_limits<double>::infinity();
      double db = std::numeric_limits<double>::infinity();

      for (size_t j = 0; j < N; ++j) {
        if (keep[ord[j]]) {
          da = std::min(da, std::abs(arc[a] - arc[j]));
          db = std::min(db, std::abs(arc[b] - arc[j]));
        }
      }

      return da > db;
    });

    for (size_t pos : pos_order) {
      size_t i = ord[pos];
      if (!keep[i]) {
        keep[i] = 1;
        ++nKeep;
        if (nKeep >= cap) break;
      }
    }
  }

  // Rebuild archive in Pareto-front order
  std::vector<arma::uvec> plans2;
  std::vector<double> U2, S2;

  plans2.reserve(cap);
  U2.reserve(cap);
  S2.reserve(cap);

  for (size_t k = 0; k < N; ++k) {
    size_t i = ord[k];
    if (keep[i]) {
      plans2.push_back(plans[i]);
      U2.push_back(U[i]);
      S2.push_back(S[i]);
    }
  }

  plans.swap(plans2);
  U.swap(U2);
  S.swap(S2);
}

static inline void inject_archive_arclength_into_parents(
    arma::umat& Parents,
    arma::vec& uParents,
    arma::vec& simParents,
    const std::vector<arma::uvec>& arch_plans,
    const std::vector<double>& arch_u,
    const std::vector<double>& arch_sim,
    arma::uword nInject,
    uint64_t seed)
{
  const size_t N = arch_u.size();
  if (N == 0 || nInject == 0) return;

  nInject = std::min<arma::uword>(nInject, Parents.n_cols);

  // Sort archive along Pareto front by similarity
  std::vector<size_t> ord(N);
  for (size_t i = 0; i < N; ++i) ord[i] = i;

  std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {
    if (arch_sim[a] != arch_sim[b]) return arch_sim[a] < arch_sim[b];
    return arch_u[a] > arch_u[b];
  });

  // Pick approximately evenly spaced archive positions
  std::vector<size_t> chosen;
  chosen.reserve(nInject);

  if (nInject == 1) {
    chosen.push_back(ord[N / 2]);
  } else {
    std::vector<double> targets = {
      0.00, 0.02, 0.05, 0.10, 0.20,
      0.80, 0.90, 0.95, 0.98, 1.00
    };

    for (arma::uword k = 0; k < nInject; ++k) {
      double target = targets[k % targets.size()];
      size_t idx = static_cast<size_t>(
        std::round(target * static_cast<double>(N - 1))
      );
      chosen.push_back(ord[idx]);
    }
  }

  // Randomize which parent slots are replaced
  std::mt19937 gen((uint32_t)seed);
  std::uniform_int_distribution<arma::uword> dis(0, Parents.n_cols - 1);

  std::vector<char> used_slot(Parents.n_cols, 0);

  for (arma::uword k = 0; k < chosen.size(); ++k) {
    arma::uword slot = dis(gen);

    // avoid replacing same parent slot twice
    for (int tries = 0; tries < 64 && used_slot[slot]; ++tries) {
      slot = dis(gen);
    }

    if (used_slot[slot]) {
      for (arma::uword s = 0; s < Parents.n_cols; ++s) {
        if (!used_slot[s]) {
          slot = s;
          break;
        }
      }
    }

    used_slot[slot] = 1;

    size_t a = chosen[k];

    Parents.col(slot) = arch_plans[a];
    uParents(slot) = arch_u[a];
    simParents(slot) = arch_sim[a];
  }
}
// ------------------------------
// Main
// ------------------------------

// [[Rcpp::export]]
Rcpp::List cpp_optimal_cross_pareto(const NumericMatrix& Crosses,
                                    const NumericMatrix& fixedCrosses,
                                    arma::uword nCross,
                                    arma::mat& G,
                                    arma::vec& u,
                                    arma::vec& ufixed,
                                    const arma::vec& minContr,
                                    const arma::vec& maxContr,
                                    const arma::mat& A,
                                    const arma::mat& Afixed,
                                    const arma::vec& lower,
                                    const arma::vec& upper,
                                    double probMut=0.01,
                                    arma::uword nMutate=0,
                                    arma::uword nSel=500,
                                    arma::uword nPop=10000,
                                    arma::uword maxGen=500,
                                    arma::uword maxRun=100,
                                    int nThreads=4)
{
  ct_set_threads(nThreads);

  const arma::uword potCross  = Crosses.nrow();
  const arma::uword fixedRows = fixedCrosses.nrow();

  if (nCross == 0) {
    Rcpp::stop("nCross must be greater than zero.");
  }

  if (fixedRows > nCross) {
    Rcpp::stop(
      "The number of fixed crosses cannot exceed nCross."
    );
  }

  const arma::uword nVar = nCross - fixedRows;

  if (nVar > potCross) {
    Rcpp::stop(
      "The number of variable crosses exceeds the number "
      "of available candidate crosses."
    );
  }


  // Validate dimensions and population settings before allocating matrices.
  if (Crosses.ncol() != 2 || fixedCrosses.ncol() != 2 ||
      G.n_rows != G.n_cols || u.n_elem != potCross ||
      ufixed.n_elem != fixedRows) {
    Rcpp::stop("Invalid Crosses, fixedCrosses, G, u, or ufixed dimensions.");
  }
  if (nSel < 2 || nPop < 2 || nSel > nPop ||
      static_cast<double>(nPop) >
          static_cast<double>(nSel) * (nSel - 1) / 2.0) {
    Rcpp::stop("Require 2 <= nSel <= nPop <= choose(nSel, 2).");
  }

  const arma::uword nInd = G.n_cols;
  for (int i = 0; i < Crosses.nrow(); ++i) {
    for (int j = 0; j < 2; ++j) {
      if (!std::isfinite(Crosses(i, j)) || Crosses(i, j) < 1 ||
          Crosses(i, j) > nInd ||
          std::floor(Crosses(i, j)) != Crosses(i, j)) {
        Rcpp::stop("Crosses must contain valid 1-based integer parent indices.");
      }
    }
  }
  for (int i = 0; i < fixedCrosses.nrow(); ++i) {
    for (int j = 0; j < 2; ++j) {
      if (!std::isfinite(fixedCrosses(i, j)) || fixedCrosses(i, j) < 1 ||
          fixedCrosses(i, j) > nInd ||
          std::floor(fixedCrosses(i, j)) != fixedCrosses(i, j)) {
        Rcpp::stop("fixedCrosses must contain valid 1-based integer parent indices.");
      }
    }
  }

  arma::umat Crosses_mat = as<arma::umat>(Crosses);
  arma::umat fixedCrosses_mat = as<arma::umat>(fixedCrosses);
  Crosses_mat -= 1;
  fixedCrosses_mat -= 1;

  // Reject individually impossible parental and cross-level constraints.
  check_constraint_feasibility(
    Crosses_mat, fixedCrosses_mat, nCross, nInd,
    minContr, maxContr, A, Afixed, lower, upper
  );

  // Contribution of fixed crosses and their numerical constraint scores.
  arma::vec xfixed = calcContr(fixedCrosses_mat, nInd, nCross);
  arma::mat scoreByConstraint = A.t();
  arma::vec fixedExtra(A.n_cols, arma::fill::zeros);
  for (arma::uword f = 0; f < Afixed.n_rows; ++f) {
    for (arma::uword j = 0; j < A.n_cols; ++j) {
      fixedExtra(j) += Afixed(f, j) / nCross;
    }
  }

  ConstraintSystem constraints{
    Crosses_mat, scoreByConstraint, minContr, maxContr,
    lower, upper, xfixed, fixedExtra, nCross, nInd, 1e-12
  };
  double uMax = 0.0, simMax = 0.0;



  // seed
  std::random_device rd;
  const uint64_t base_seed = ((uint64_t)rd() << 32) ^ (uint64_t)rd();

  const double ufixedSum = arma::accu(ufixed);

  // GA containers
  arma::umat Progeny(nVar, nPop, arma::fill::zeros);
  arma::umat Parents(nVar, nSel, arma::fill::zeros);
  arma::uvec Best(nVar, arma::fill::zeros);

  arma::vec uProgeny(nPop, arma::fill::zeros), uParents(nSel, arma::fill::zeros);
  arma::vec simProgeny(nPop, arma::fill::zeros), simParents(nSel, arma::fill::zeros);

  arma::uvec rankProgeny(nPop, arma::fill::zeros);

  // --- baseline uMax/simMax ---
  arma::uvec uMax_plan;

  // Keep the complete ranking because the constrained procedure
  // may need crosses below the first nVar positions.
  arma::uvec uOrder =
    arma::sort_index(u, "descend");

  // High-gain anchor must satisfy parental and all cross-level constraints.
  uMax_plan = build_constrained_umax_plan(
    u, uOrder, constraints, nVar, base_seed
  );

  {
    arma::vec x =
      calcContr(
        Crosses_mat.rows(uMax_plan),
        nInd,
        nCross
      );

    uMax =
      (arma::accu(u.elem(uMax_plan)) + ufixedSum) /
        static_cast<double>(nCross);

    simMax =
      calc_sim(x + xfixed, G);
  }


  // =========================================================
  // Phase 1: GA minimizing sim only
  // =========================================================
  Rcpp::Rcout << "\nPhase 1: Optimize for Similarity (minimize sim)\n";

  // Generate initial plans
  ct_parallel_for(0, (int)nPop, [&](int i) {
    std::mt19937 gen(
        static_cast<uint32_t>(
          base_seed ^
            (0x9E3779B97F4A7C15ULL +
            static_cast<uint64_t>(i))
        )
    );

    Progeny.col((arma::uword)i) =
      sampleInt_std(nVar, potCross, gen);
  });

  // Repair all initial plans
  repair_constraint_plans(
    Progeny, constraints, uMax_plan,
    base_seed ^ 0x1010101010101010ULL
  );

  // Evaluate the repaired plans
  ct_parallel_for(0, (int)nPop, [&](int i) {
    arma::uvec idx =
      arma::conv_to<arma::uvec>::from(
        Progeny.col((arma::uword)i)
      );

    double uu, ss;

    eval_plan(
      idx,
      u,
      ufixedSum,
      nCross,
      G,
      Crosses_mat,
      nInd,
      xfixed,
      uu,
      ss
    );

    simProgeny((arma::uword)i) = ss;
  });


  rankProgeny = arma::sort_index(simProgeny, "ascend");
  for (arma::uword i = 0; i < nSel; ++i) {
    Parents.col(i) = Progeny.col(rankProgeny(i));
    simParents(i)  = simProgeny(rankProgeny(i));
  }

  double simBest_phase1 = simParents(0);
  arma::uvec Best_sim_phase1 = arma::conv_to<arma::uvec>::from(Parents.col(0)); // copy
  arma::uword stall_sim = 0;

  Rcpp::Rcout << "Gen  simBest\n";
  for (arma::uword gen = 0; gen < maxGen; ++gen) {
    std::mt19937 rng_plan((uint32_t)(base_seed ^ (0xC3A5C85C97CB3127ULL + (uint64_t)gen)));
    arma::umat crossPlan = sampHalfDialComb_std(nSel, nPop, rng_plan);

    // Generate children only
    ct_parallel_for(0, (int)nPop, [&](int i) {
      std::mt19937 gen_i(
          static_cast<uint32_t>(
            base_seed ^
              (0xD2B74407B1CE6E93ULL +
              static_cast<uint64_t>(gen) * 1315423911ULL +
              static_cast<uint64_t>(i))
          )
      );

      arma::uvec child = mate_uniform(
        arma::conv_to<arma::uvec>::from(
          Parents.col(
            crossPlan(0, (arma::uword)i)
          )
        ),
        arma::conv_to<arma::uvec>::from(
          Parents.col(
            crossPlan(1, (arma::uword)i)
          )
        ),
        gen_i,
        potCross
      );

      std::uniform_real_distribution<> U01(0.0, 1.0);

      if (U01(gen_i) < probMut) {
        child =
          mutate(
            child,
            nMutate,
            potCross,
            gen_i
          );
      }

      Progeny.col((arma::uword)i) = child;
    });

    // Repair complete population once
    repair_constraint_plans(
    Progeny, constraints, uMax_plan,
    base_seed ^
        (0x2020202020202020ULL +
        static_cast<uint64_t>(gen))
  );

    // Evaluate repaired population
    ct_parallel_for(0, (int)nPop, [&](int i) {
      arma::uvec idx =
        arma::conv_to<arma::uvec>::from(
          Progeny.col((arma::uword)i)
        );

      double uu, ss;

      eval_plan(
        idx,
        u,
        ufixedSum,
        nCross,
        G,
        Crosses_mat,
        nInd,
        xfixed,
        uu,
        ss
      );

      simProgeny((arma::uword)i) = ss;
    });


    rankProgeny = arma::sort_index(simProgeny, "ascend");
    for (arma::uword i = 0; i < nSel; ++i) {
      Parents.col(i) = Progeny.col(rankProgeny(i));
      simParents(i)  = simProgeny(rankProgeny(i));
    }

    if (simParents(0) < simBest_phase1) {
      simBest_phase1 = simParents(0);
      Best_sim_phase1 = arma::conv_to<arma::uvec>::from(Parents.col(0)); // copy
      stall_sim = 0;
    } else {
      ++stall_sim;
    }

    if (gen % 10 == 0) Rcpp::Rcout << gen << "  " << simBest_phase1 << "\n";
    if (stall_sim >= maxRun) break;
  }

  // simPool = best low-sim from final Progeny
  arma::uvec ordSim = arma::sort_index(simProgeny, "ascend");
  arma::umat simPool(nVar, nPop);
  for (arma::uword i = 0; i < nPop; ++i) simPool.col(i) = Progeny.col(ordSim(i));

  // =========================================================
  // Init: epsilon-constraint sweep (maximize u s.t. sim <= eps)
  // =========================================================
  Rcpp::Rcout << "\nInit: epsilon-constraint sweep (maximize u s.t. sim<=eps)\n";

  // equally spaced eps from best(sim-only) -> simMax
  const arma::uword N_EPS = 20; // tune
  double eps_lo = simBest_phase1;
  double eps_hi = simMax;
  if (!(eps_hi >= eps_lo)) eps_hi = eps_lo;

  std::vector<double> eps_grid;
  eps_grid.reserve((size_t)N_EPS);
  if (N_EPS <= 1) {
    eps_grid.push_back(eps_lo);
  } else {
    for (arma::uword k = 0; k < N_EPS; ++k) {
      double t = (double)k / (double)(N_EPS - 1);
      eps_grid.push_back(eps_lo + t * (eps_hi - eps_lo));
    }
  }

  auto better_eps = [&](double uA, double sA, double eps,
                        double uB, double sB) -> bool {
                          const double vA = std::max(0.0, sA - eps);
                          const double vB = std::max(0.0, sB - eps);
                          const bool fA = (vA <= 0.0);
                          const bool fB = (vB <= 0.0);

                          if (fA != fB) return fA;
                          if (fA && fB) {
                            if (uA != uB) return uA > uB;
                            return sA < sB;
                          }
                          if (vA != vB) return vA < vB;
                          if (uA != uB) return uA > uB;
                          return sA < sB;
                        };

  // store best plan per eps (copies)
  std::vector<arma::uvec> eps_best_plans;
  eps_best_plans.reserve(eps_grid.size() + 8);

  // anchors
  if (Best_sim_phase1.n_elem == nVar) eps_best_plans.push_back(Best_sim_phase1);
  if (uMax_plan.n_elem == nVar) eps_best_plans.push_back(uMax_plan);

  for (size_t wi = 0; wi < eps_grid.size(); ++wi) {
    const double eps = eps_grid[wi];

    Rcpp::Rcout << "\n  Eps run " << (wi + 1) << "/" << eps_grid.size()
                << "  eps=" << eps << "\n";

    // init random
    ct_parallel_for(0, (int)nPop, [&](int i) {
      std::mt19937 gen((uint32_t)(base_seed ^ (0xA5A5A5A55A5A5A5AULL
                                                  + (uint64_t)i
                                                  + (uint64_t)wi * 0x9E3779B1ULL)));
                                                  Progeny.col((arma::uword)i) = sampleInt_std(nVar, potCross, gen);
    });

    // inject anchors + some user seeds
    arma::uword inj = 0;
    if (inj < nPop && Best_sim_phase1.n_elem == nVar) Progeny.col(inj++) = Best_sim_phase1;
    if (inj < nPop && uMax_plan.n_elem == nVar)       Progeny.col(inj++) = uMax_plan;




    repair_constraint_plans(
    Progeny, constraints, uMax_plan,
    base_seed ^
        (0x3030303030303030ULL +
        static_cast<uint64_t>(wi))
  );
    enforce_unique_population(Progeny, constraints, uMax_plan, base_seed);

    // eval init (u/sim only)
    ct_parallel_for(0, (int)nPop, [&](int i) {
      arma::uvec idx = arma::conv_to<arma::uvec>::from(Progeny.col((arma::uword)i));
      double uu, ss;
      eval_plan(idx, u, ufixedSum, nCross, G, Crosses_mat, nInd, xfixed, uu, ss);

      uProgeny((arma::uword)i)   = uu;
      simProgeny((arma::uword)i) = ss;
    });

    // select parents by feasibility-first
    std::vector<int> ord((size_t)nPop);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b){
      return better_eps(uProgeny((arma::uword)a), simProgeny((arma::uword)a), eps,
                        uProgeny((arma::uword)b), simProgeny((arma::uword)b));
    });

    for (arma::uword i = 0; i < nSel; ++i) {
      Parents.col(i) = Progeny.col((arma::uword)ord[(size_t)i]);
      uParents(i)    = uProgeny((arma::uword)ord[(size_t)i]);
      simParents(i)  = simProgeny((arma::uword)ord[(size_t)i]);
    }

    double best_u   = uParents(0);
    double best_sim = simParents(0);
    double best_viol = std::max(0.0, best_sim - eps);
    int best_feas = (best_viol <= 0.0) ? 1 : 0;
    arma::uvec best_plan_eps = arma::conv_to<arma::uvec>::from(Parents.col(0)); // copy

    arma::uword stall_eps = 0;

    Rcpp::Rcout << "  Gen  best_u  best_sim  feas  viol\n";
    Rcpp::Rcout << "  " << 0 << "  " << best_u << "  " << best_sim
                << "  " << best_feas << "  " << best_viol << "\n";

    for (arma::uword gen = 0; gen < maxGen; ++gen) {
      std::mt19937 rng_plan((uint32_t)(base_seed ^ (0xC3A5C85C97CB3127ULL
                                                       + (uint64_t)gen
                                                       + (uint64_t)wi * 0xD1B54A32ULL)));
                                                       arma::umat crossPlan = sampHalfDialComb_std(nSel, nPop, rng_plan);

                                                       // Generate children
                                                       ct_parallel_for(0, (int)nPop, [&](int i) {
                                                         std::mt19937 gen_i(
                                                             static_cast<uint32_t>(
                                                               base_seed ^
                                                                 (0x9E3779B97F4A7C15ULL +
                                                                 static_cast<uint64_t>(gen) * 1315423911ULL +
                                                                 static_cast<uint64_t>(i) +
                                                                 static_cast<uint64_t>(wi) * 0x85EBCA6BULL)
                                                             )
                                                         );

                                                         arma::uvec child = mate_uniform(
                                                           arma::conv_to<arma::uvec>::from(
                                                             Parents.col(
                                                               crossPlan(0, (arma::uword)i)
                                                             )
                                                           ),
                                                           arma::conv_to<arma::uvec>::from(
                                                             Parents.col(
                                                               crossPlan(1, (arma::uword)i)
                                                             )
                                                           ),
                                                           gen_i,
                                                           potCross
                                                         );

                                                         std::uniform_real_distribution<> U01(0.0, 1.0);

                                                         if (U01(gen_i) < probMut) {
                                                           child =
                                                             mutate(
                                                               child,
                                                               nMutate,
                                                               potCross,
                                                               gen_i
                                                             );
                                                         }

                                                         Progeny.col((arma::uword)i) = child;
                                                       });

                                                       // Repair children
                                                       repair_constraint_plans(
    Progeny, constraints, uMax_plan,
    base_seed ^
                                                           (0x4040404040404040ULL +
                                                           static_cast<uint64_t>(wi) * 1000003ULL +
                                                           static_cast<uint64_t>(gen))
  );

                                                       // Evaluate repaired children
                                                       ct_parallel_for(0, (int)nPop, [&](int i) {
                                                         arma::uvec idx =
                                                           arma::conv_to<arma::uvec>::from(
                                                             Progeny.col((arma::uword)i)
                                                           );

                                                         double uu, ss;

                                                         eval_plan(
                                                           idx,
                                                           u,
                                                           ufixedSum,
                                                           nCross,
                                                           G,
                                                           Crosses_mat,
                                                           nInd,
                                                           xfixed,
                                                           uu,
                                                           ss
                                                         );

                                                         uProgeny((arma::uword)i)   = uu;
                                                         simProgeny((arma::uword)i) = ss;
                                                       });

                                                       std::iota(ord.begin(), ord.end(), 0);
                                                       std::sort(ord.begin(), ord.end(), [&](int a, int b){
                                                         return better_eps(uProgeny((arma::uword)a), simProgeny((arma::uword)a), eps,
                                                                           uProgeny((arma::uword)b), simProgeny((arma::uword)b));
                                                       });

                                                       for (arma::uword i = 0; i < nSel; ++i) {
                                                         Parents.col(i) = Progeny.col((arma::uword)ord[(size_t)i]);
                                                         uParents(i)    = uProgeny((arma::uword)ord[(size_t)i]);
                                                         simParents(i)  = simProgeny((arma::uword)ord[(size_t)i]);
                                                       }

                                                       bool improved = better_eps(uParents(0), simParents(0), eps, best_u, best_sim);
                                                       if (improved) {
                                                         best_u   = uParents(0);
                                                         best_sim = simParents(0);
                                                         best_viol = std::max(0.0, best_sim - eps);
                                                         best_feas = (best_viol <= 0.0) ? 1 : 0;
                                                         best_plan_eps = arma::conv_to<arma::uvec>::from(Parents.col(0)); // copy
                                                         stall_eps = 0;
                                                       } else {
                                                         ++stall_eps;
                                                       }

                                                       if (gen % 10 == 0) {
                                                         Rcpp::Rcout << "  " << (gen + 1) << "  " << best_u << "  " << best_sim
                                                                     << "  " << best_feas << "  " << best_viol << "\n";
                                                       }
                                                       if (stall_eps >= maxRun) break;
    }

    if (best_plan_eps.n_elem == nVar) eps_best_plans.push_back(best_plan_eps);
  }

  // make epsPool unique (by sorted key)
  std::unordered_set<std::string> eps_seen;
  eps_seen.reserve(eps_best_plans.size() * 2);

  std::vector<arma::uvec> eps_uniq;
  eps_uniq.reserve(eps_best_plans.size());

  for (size_t i = 0; i < eps_best_plans.size(); ++i) {
    if (eps_best_plans[i].n_elem != nVar) continue;
    arma::uvec vs = arma::sort(eps_best_plans[i]);
    std::string key = plan_key_sorted(vs);
    if (eps_seen.insert(key).second) eps_uniq.push_back(eps_best_plans[i]);
  }

  arma::uword nEpsPool = (arma::uword)eps_uniq.size();
  arma::umat epsPool(nVar, nEpsPool, arma::fill::zeros);
  for (arma::uword i = 0; i < nEpsPool; ++i) epsPool.col(i) = eps_uniq[(size_t)i];

  // =========================================================
  // Init: build fused pool then Pareto-filter only
  // Sources: simPool, epsPool, random, user seeds
  // =========================================================
  Rcpp::Rcout << "\nInit: build fused pool then Pareto-filter only\n";

  const arma::uword nSimUse  = std::min<arma::uword>(nPop, simPool.n_cols);
  const arma::uword nEpsUse  = std::min<arma::uword>(nPop, epsPool.n_cols); // cap to avoid huge init
  const arma::uword nRanUse  = nPop;

  const arma::uword nBig = nSimUse + nEpsUse + nRanUse;

  arma::umat ProgenyBig(nVar, nBig, arma::fill::zeros);

  // debug snapshots
  arma::umat simInit(nVar, nBig, arma::fill::zeros);
  arma::umat epsInit(nVar, nBig, arma::fill::zeros);
  arma::umat ranInit(nVar, nBig, arma::fill::zeros);

  arma::uword used = 0;

  // 1) simPool
  for (arma::uword k = 0; k < nSimUse; ++k) {
    ProgenyBig.col(used) = simPool.col(k);
    simInit.col(used)    = ProgenyBig.col(used);
    ++used;
  }

  // 2) epsPool
  for (arma::uword k = 0; k < nEpsUse; ++k) {
    ProgenyBig.col(used) = epsPool.col(k);
    epsInit.col(used)    = ProgenyBig.col(used);
    ++used;
  }

  // 3) random
  ct_parallel_for((int)used, (int)(used + nRanUse), [&](int i) {
    std::mt19937 gen_i((uint32_t)(base_seed ^ (0xB492B66FBE98F273ULL + (uint64_t)i)));
    ProgenyBig.col((arma::uword)i) = sampleInt_std(nVar, potCross, gen_i);
  });
  for (arma::uword i = used; i < used + nRanUse; ++i) ranInit.col(i) = ProgenyBig.col(i);
  used += nRanUse;



  if (used != nBig) {
    Rcpp::stop("Internal error: used != nBig in fused initialization.");
  }



  repair_constraint_plans(
    ProgenyBig, constraints, uMax_plan,
    base_seed ^ 0x5050505050505050ULL
  );
  enforce_unique_population(ProgenyBig, constraints, uMax_plan, base_seed);

  // Evaluate big pool with the REAL objective weights (weightu/weightcrossover)
  arma::vec uBig(nBig, arma::fill::zeros);
  arma::vec simBig(nBig, arma::fill::zeros);

  ct_parallel_for(0, (int)nBig, [&](int i) {
    arma::uvec idx = arma::conv_to<arma::uvec>::from(ProgenyBig.col((arma::uword)i));
    double uu, ss;

    eval_plan(idx, u, ufixedSum, nCross, G, Crosses_mat, nInd, xfixed, uu, ss);

    uBig((arma::uword)i)   = uu;
    simBig((arma::uword)i) = ss;
  });

  // Pareto-filter indices among big pool (maximize u, minimize sim)
  std::vector<int> big_pareto;
  big_pareto.reserve((size_t)nBig);

  {
    std::vector<char> keep((size_t)nBig, 1);
    for (arma::uword i = 0; i < nBig; ++i) {
      if (!keep[(size_t)i]) continue;
      for (arma::uword j = 0; j < nBig; ++j) {
        if (i == j || !keep[(size_t)j]) continue;
        if (dominates_u_sim(uBig(j), simBig(j), uBig(i), simBig(i))) { keep[(size_t)i] = 0; break; }
      }
    }
    for (arma::uword i = 0; i < nBig; ++i) if (keep[(size_t)i]) big_pareto.push_back((int)i);
  }

  // choose init population indices: ONLY from Pareto if possible
  std::vector<int> init_idx;
  init_idx.reserve((size_t)nPop);

  if ((int)big_pareto.size() >= (int)nPop) {
    // downselect *within Pareto only* (crowding diversity)
    init_idx = nsga2_select_maxu_minsim(uBig, simBig, big_pareto, (int)nPop);
  } else {
    // take all Pareto, then fill remainder by best u (tie lower sim)
    init_idx = big_pareto;

    std::vector<int> rest;
    rest.reserve((size_t)nBig - init_idx.size());
    std::vector<char> in_pf((size_t)nBig, 0);
    for (int x : big_pareto) in_pf[(size_t)x] = 1;
    for (arma::uword i = 0; i < nBig; ++i) if (!in_pf[(size_t)i]) rest.push_back((int)i);

    std::sort(rest.begin(), rest.end(), [&](int a, int b){
      if (uBig((arma::uword)a) != uBig((arma::uword)b)) return uBig((arma::uword)a) > uBig((arma::uword)b);
      return simBig((arma::uword)a) < simBig((arma::uword)b);
    });

    for (size_t k = 0; k < rest.size() && init_idx.size() < (size_t)nPop; ++k) {
      init_idx.push_back(rest[k]);
    }
  }

  // Write selected to Progeny (nVar x nPop)
  for (arma::uword i = 0; i < nPop; ++i) {
    int j = init_idx[(size_t)i];
    arma::uvec plan = arma::conv_to<arma::uvec>::from(ProgenyBig.col((arma::uword)j));

    Progeny.col(i) = plan;

    uProgeny(i)   = uBig((arma::uword)j);
    simProgeny(i) = simBig((arma::uword)j);
  }

  // uniqueness might change plans -> re-evaluate


  repair_constraint_plans(
    Progeny, constraints, uMax_plan,
    base_seed ^ 0x6060606060606060ULL
  );
  enforce_unique_population(Progeny, constraints, uMax_plan, base_seed);

  ct_parallel_for(0, (int)nPop, [&](int i) {
    arma::uvec idx = arma::conv_to<arma::uvec>::from(Progeny.col((arma::uword)i));
    double uu, ss;
    eval_plan(idx, u, ufixedSum, nCross, G, Crosses_mat, nInd, xfixed, uu, ss);

    uProgeny((arma::uword)i)   = uu;
    simProgeny((arma::uword)i) = ss;
  });

  // =========================================================
  // Archive + main NSGA-II GA loop (same as your remainder)
  // =========================================================
  std::vector<arma::uvec> arch_plans;
  std::vector<double> arch_u, arch_sim;
  const size_t ARCH_CAP = 5000;

  arch_plans.reserve((size_t)nPop * 4);
  arch_u.reserve((size_t)nPop * 4);
  arch_sim.reserve((size_t)nPop * 4);

  for (arma::uword i = 0; i < nPop; ++i) {
    arch_plans.push_back(arma::conv_to<arma::uvec>::from(Progeny.col(i)));
    arch_u.push_back(uProgeny(i));
    arch_sim.push_back(simProgeny(i));
  }
  pareto_filter_archive(arch_plans, arch_u, arch_sim);
  if (ARCH_CAP > 0) cap_archive_by_arc_length(arch_plans, arch_u, arch_sim, ARCH_CAP);

  // select parents via NSGA-II on FULL population
  std::vector<int> cand((size_t)nPop);
  std::iota(cand.begin(), cand.end(), 0);
  std::vector<int> sel = nsga2_select_maxu_minsim(uProgeny, simProgeny, cand, (int)nSel);

  for (arma::uword i = 0; i < nSel; ++i) {
    int idx = sel[(size_t)i];
    Parents.col(i) = Progeny.col((arma::uword)idx);
    uParents(i)    = uProgeny((arma::uword)idx);
    simParents(i)  = simProgeny((arma::uword)idx);
  }

  // pick best for logging: highest u, tie-break lowest sim
  arma::uword best_i = 0;
  for (arma::uword i = 1; i < nSel; ++i) {
    if (uParents(i) > uParents(best_i) ||
        (uParents(i) == uParents(best_i) && simParents(i) < simParents(best_i))) {
      best_i = i;
    }
  }

  Best = arma::conv_to<arma::uvec>::from(Parents.col(best_i));
  double uBest_val = uParents(best_i);
  double simBest   = simParents(best_i);

  // GA trace
  std::vector<double> trace_gen, trace_best_u, trace_best_sim;
  std::vector<double> trace_mean_u, trace_mean_sim;

  auto log_generation = [&](arma::uword gen_index) {
    trace_gen.push_back((double)gen_index);
    trace_best_u.push_back(uBest_val);
    trace_best_sim.push_back(simBest);
    trace_mean_u.push_back(arma::mean(uProgeny));
    trace_mean_sim.push_back(arma::mean(simProgeny));
  };

  log_generation(0);

  Rcpp::Rcout << "Gen  uBest  simBest\n";
  Rcpp::Rcout << 0 << "  " << uBest_val << "  " << simBest << "\n";


  arma::uword stall_main = 0;

  for (arma::uword gen = 0; gen < maxGen; ++gen) {
    std::mt19937 gen_plan((uint32_t)(base_seed ^ (0x9E3779B185EBCA87ULL + (uint64_t)gen)));
    arma::umat crossPlan = sampHalfDialComb_std(nSel, nPop, gen_plan);

    // Generate children only
    ct_parallel_for(0, (int)nPop, [&](int i) {
      std::mt19937 gen_i(
          static_cast<uint32_t>(
            base_seed ^
              (0xC949D7C7509E6557ULL +
              static_cast<uint64_t>(gen) * 40503ULL +
              static_cast<uint64_t>(i))
          )
      );

      arma::uvec child = mate_uniform(
        arma::conv_to<arma::uvec>::from(
          Parents.col(
            crossPlan(0, (arma::uword)i)
          )
        ),
        arma::conv_to<arma::uvec>::from(
          Parents.col(
            crossPlan(1, (arma::uword)i)
          )
        ),
        gen_i,
        potCross
      );

      std::uniform_real_distribution<> U01(0.0, 1.0);

      if (U01(gen_i) < probMut) {
        child =
          mutate(
            child,
            nMutate,
            potCross,
            gen_i
          );
      }

      Progeny.col((arma::uword)i) = child;
    });

    // Repair population exactly once
    repair_constraint_plans(
    Progeny, constraints, uMax_plan,
    base_seed ^
        (0x7070707070707070ULL +
        static_cast<uint64_t>(gen))
  );

    // Evaluate repaired plans
    ct_parallel_for(0, (int)nPop, [&](int i) {
      arma::uvec idx =
        arma::conv_to<arma::uvec>::from(
          Progeny.col((arma::uword)i)
        );

      double uu, ss;

      eval_plan(
        idx,
        u,
        ufixedSum,
        nCross,
        G,
        Crosses_mat,
        nInd,
        xfixed,
        uu,
        ss
      );

      uProgeny((arma::uword)i)   = uu;
      simProgeny((arma::uword)i) = ss;
    });

    // archive update
    for (arma::uword i = 0; i < nPop; ++i) {
      arch_plans.push_back(arma::conv_to<arma::uvec>::from(Progeny.col(i)));
      arch_u.push_back(uProgeny(i));
      arch_sim.push_back(simProgeny(i));
    }
    pareto_filter_archive(arch_plans, arch_u, arch_sim);
    if (ARCH_CAP > 0) cap_archive_by_arc_length(arch_plans, arch_u, arch_sim, ARCH_CAP);

    // parent selection
    std::vector<int> cand2((size_t)nPop);
    std::iota(cand2.begin(), cand2.end(), 0);
    std::vector<int> sel2 = nsga2_select_maxu_minsim(uProgeny, simProgeny, cand2, (int)nSel);

    for (arma::uword i = 0; i < nSel; ++i) {
        const int idx = sel2[(size_t)i];

        Parents.col(i) = Progeny.col((arma::uword)idx);
        uParents(i)    = uProgeny((arma::uword)idx);
        simParents(i)  = simProgeny((arma::uword)idx);
      }

      const arma::uword nInject =
        std::max<arma::uword>(
          2,
          (nSel / 100) * 25
        );

      if (gen % 5 == 0) {
        inject_archive_arclength_into_parents(
          Parents,
          uParents,
          simParents,
          arch_plans,
          arch_u,
          arch_sim,
          nInject,
          base_seed ^
            (0x123456789ABCDEFULL +
            static_cast<uint64_t>(gen))
        );
      }

    best_i = 0;
    for (arma::uword i = 1; i < nSel; ++i) {
      if (uParents(i) > uParents(best_i) ||
          (uParents(i) == uParents(best_i) && simParents(i) < simParents(best_i))) {
        best_i = i;
      }
    }

    // keep previous best (u, sim) for stall detection
    double prev_u = uBest_val;
    double prev_s = simBest;

    // update current best from selected parents
    uBest_val = uParents(best_i);
    simBest   = simParents(best_i);
    Best      = arma::conv_to<arma::uvec>::from(Parents.col(best_i));

    // stall logic: reset if improved in u, or same u but better (smaller) sim
    if (uBest_val > prev_u || (uBest_val == prev_u && simBest < prev_s)) {
      stall_main = 0;
    } else {
      ++stall_main;
    }

    log_generation(gen + 1);

    if (gen % 10 == 0) {
      Rcpp::Rcout << (gen + 1) << "  " << uBest_val << "  " << simBest << "\n";
    }
    if (stall_main >= maxRun) break;
  }

  // Final best plan
  Best = arma::sort(Best);
  const ConstraintState finalState = constraints.initialize(Best);
  if (finalState.violation > 1e-10) {
    Rcpp::stop("Internal error: best Pareto plan violates constraints.");
  }
  arma::umat outCrossPlan = arma::join_cols(Crosses_mat.rows(Best), fixedCrosses_mat);

  // Final Pareto outputs from archive
  Rcpp::DataFrame pareto_df;
  Rcpp::List pareto_crossplans;
  Rcpp::NumericMatrix pareto_constraint_values;
  {
    const size_t nPF = arch_u.size();

    std::vector<size_t> ord(nPF);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b){
      return arch_sim[a] < arch_sim[b];
    });
    Rcpp::Rcout << "optimization finished \n";

    pareto_crossplans = Rcpp::List((int)nPF);
    Rcpp::NumericVector pu((int)nPF), ps((int)nPF);
    pareto_constraint_values = Rcpp::NumericMatrix(
      static_cast<int>(nPF), static_cast<int>(A.n_cols));

    for (size_t k = 0; k < nPF; ++k) {
      size_t a = ord[k];
      pu[(int)k] = arch_u[a];
      ps[(int)k] = arch_sim[a];

      arma::uvec vv = arma::sort(arch_plans[a]);
      const ConstraintState st = constraints.initialize(vv);
      if (st.violation > 1e-10) {
        Rcpp::stop("Internal error: infeasible plan found in Pareto archive.");
      }
      for (arma::uword j = 0; j < A.n_cols; ++j) {
        pareto_constraint_values(static_cast<int>(k), static_cast<int>(j)) = st.extra(j);
      }
      arma::umat cp = arma::join_cols(Crosses_mat.rows(vv), fixedCrosses_mat);
      pareto_crossplans[(int)k] = cp + 1;
    }

    pareto_df = Rcpp::DataFrame::create(
      _["pareto_id"] = Rcpp::seq(1, (int)nPF),
      _["u"]   = pu,
      _["sim"] = ps
    );
  }

  // GA trace to DataFrame
  int nTrace = (int)trace_gen.size();
  NumericVector genVec(nTrace), best_u(nTrace), best_sim_v(nTrace),
  mean_u(nTrace), mean_sim_v(nTrace);

  for (int i = 0; i < nTrace; ++i) {
    genVec[i]       = trace_gen[i];
    best_u[i]       = trace_best_u[i];
    best_sim_v[i]   = trace_best_sim[i];
    mean_u[i]       = trace_mean_u[i];
    mean_sim_v[i]   = trace_mean_sim[i];
  }

  DataFrame gaTrace = DataFrame::create(
    _["generation"]       = genVec,
    _["best_u"]           = best_u,
    _["best_sim"]         = best_sim_v,
    _["mean_u"]           = mean_u,
    _["mean_sim"]         = mean_sim_v
  );

  return Rcpp::List::create(
    Rcpp::Named("crossPlan")    = outCrossPlan + 1,
    Rcpp::Named("uMax")         = uMax,
    Rcpp::Named("simMax")       = simMax,
    Rcpp::Named("uBest")        = uBest_val,
    Rcpp::Named("simBest")      = simBest,
    Rcpp::Named("gaTrace")      = gaTrace,
    Rcpp::Named("pareto")       = pareto_df,
    Rcpp::Named("paretoPlans")  = pareto_crossplans,
    Rcpp::Named("paretoConstraintValues") = pareto_constraint_values,
    Rcpp::Named("parentContributions") = finalState.parents,
    Rcpp::Named("constraintValues") = finalState.extra
  );
}
