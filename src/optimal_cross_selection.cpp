// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <numeric>
#include <sstream>
#include <limits>
#include "parallel_backend.h"
using namespace Rcpp;


// ------------------------------
// Helpers
// ------------------------------

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



// Map linear index -> row/col in upper triangle of n x n (row <= col)
static arma::uword mapRow(const arma::uword& k, const arma::uword& n) {
  return n - 2 - static_cast<arma::uword>(
      std::floor(std::sqrt(-8.0 * static_cast<double>(k)
                             + 4.0 * static_cast<double>(n) * (static_cast<double>(n) - 1.0) - 7.0) / 2.0 - 0.5));
}
static arma::uword mapCol(const arma::uword& row, const arma::uword& k, const arma::uword& n) {
  return k + row + 1 - n * (n - 1) / 2 + (n - row) * ((n - row) - 1) / 2;
}

// Thread-safe sampling without replacement: returns sorted unique indices in [0, N-1]
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
// Half-diallel sample (n pairs from all combinations, without replacement)
static arma::umat sampHalfDialComb_std(arma::uword nLevel, arma::uword n, std::mt19937& gen) {
  const arma::uword N = nLevel * (nLevel - 1) / 2; // total possible pairs
  arma::uvec samples = sampleInt_std(std::min(n, N), N, gen); // 0..N-1
  arma::umat out(2, samples.n_elem);
  for (arma::uword i = 0; i < samples.n_elem; ++i) {
    arma::uword r = mapRow(samples(i), nLevel);
    arma::uword c = mapCol(r, samples(i), nLevel);
    out(0, i) = r;
    out(1, i) = c;
  }
  return out;
}

// Parental contributions for a crossing plan (indices are 0-based)
static arma::vec calcContr(
    const arma::umat& crosses,
    const arma::uword& nInd,
    const arma::uword& nCross)
{
  const double val =
    1.0 / (2.0 * static_cast<double>(nCross));

  arma::vec x(nInd, arma::fill::zeros);

  for (arma::uword r = 0; r < crosses.n_rows; ++r) {
    const arma::uword a = crosses(r, 0);
    const arma::uword b = crosses(r, 1);

    x(a) += val;
    x(b) += val;
  }

  return x;
}



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


// Normalize to unit square, compute angle and length
// toward high gain and low similarity.
static void calcVec(
    double& angle,
    double& length,
    double u,
    double sim,
    const double& uMax,
    const double& simMax,
    const double& uMin,
    const double& simMin
) {
  
  // Avoid division by zero when the reference values coincide.
  u = std::abs(uMax - uMin) > 1e-12
  ? (u - uMin) / (uMax - uMin)
    : 0.0;
  
  sim = std::abs(simMax - simMin) > 1e-12
  ? (simMax - sim) / (simMax - simMin)
    : 0.0;
  
  length = std::sqrt(u * u + sim * sim);
  
  angle = std::acos(
    u / std::max(length, 1e-16)
  );
  
  if (u < 0) {
    length = -length;
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
// ------------------------------
// Genetic Algorithm main
// ------------------------------

// [[Rcpp::export]]
Rcpp::List cpp_optimal_cross_selection(
    const NumericMatrix& Crosses,
    const NumericMatrix& fixedCrosses,
    arma::uword nCross,
    double targetAngle,
    arma::vec& u,
    arma::vec& ufixed,
    arma::mat& G,
    const arma::vec& minContr,
    const arma::vec& maxContr,
    const arma::mat& A,
    const arma::mat& Afixed,
    const arma::vec& lower,
    const arma::vec& upper,
    double probMut=0.01,
    arma::uword nMutate=2,
    arma::uword nSel=500,
    arma::uword nPop=10000,
    arma::uword maxGen=1000,
    arma::uword maxRun=100,
    double anglePenalty=0.5,
    int nThreads=4){
  //  threading
   
  ct_set_threads(nThreads);
  
  
  // Dimensions and input validation
  
  const arma::uword potCross =
    Crosses.nrow();
  
  const arma::uword fixedRows =
    fixedCrosses.nrow();
  
  if (
      nCross == 0 ||
        fixedRows > nCross ||
        nCross - fixedRows > potCross
  ) {
    Rcpp::stop(
      "Invalid nCross/fixedRows/candidate counts."
    );
  }
  
  const arma::uword nVar =
    nCross - fixedRows;
  
  if (
      Crosses.ncol() != 2 ||
        fixedCrosses.ncol() != 2 ||
        G.n_rows != G.n_cols ||
        u.n_elem != potCross ||
        ufixed.n_elem != fixedRows
  ) {
    Rcpp::stop(
      "Invalid crosses, G, u, or ufixed dimensions."
    );
  }
  
  if (
      nSel < 2 ||
        nPop < 2 ||
        nSel > nPop ||
        static_cast<double>(nPop) >
        static_cast<double>(nSel) * (nSel - 1) / 2.0
  ) {
    Rcpp::stop(
      "Require 2 <= nSel <= nPop <= choose(nSel, 2)."
    );
  }
  
  // Validate 1-based parent indices before converting to unsigned integers.
  for (int i = 0; i < Crosses.nrow(); ++i) {
    for (int j = 0; j < 2; ++j) {
      
      if (
          !std::isfinite(Crosses(i, j)) ||
            Crosses(i, j) < 1 ||
            Crosses(i, j) > G.n_cols ||
            std::floor(Crosses(i, j)) != Crosses(i, j)
      ) {
        Rcpp::stop(
          "Crosses must contain valid 1-based integer parent indices."
        );
      }
    }
  }
  
  
  for (int i = 0; i < fixedCrosses.nrow(); ++i) {
    for (int j = 0; j < 2; ++j) {
      
      if (
          !std::isfinite(fixedCrosses(i, j)) ||
            fixedCrosses(i, j) < 1 ||
            fixedCrosses(i, j) > G.n_cols ||
            std::floor(fixedCrosses(i, j)) != fixedCrosses(i, j)
      ) {
        Rcpp::stop(
          "fixedCrosses must contain valid 1-based integer parent indices."
        );
      }
    }
  }
  
  
  // Convert parent indices to zero-based indexing

  arma::umat Crosses_mat =
    as<arma::umat>(Crosses);
  
  arma::umat fixedCrosses_mat =
    as<arma::umat>(fixedCrosses);
  
  Crosses_mat -= 1;
  fixedCrosses_mat -= 1;
  
  arma::uword nInd =
    G.n_cols;
  
  
  // Individual feasibility checks
  // Run before allocating large GA populations.
  
  check_constraint_feasibility(
    Crosses_mat,
    fixedCrosses_mat,
    nCross,
    nInd,
    minContr,
    maxContr,
    A,
    Afixed,
    lower,
    upper
  );

  arma::umat Progeny(
      nVar,
      nPop,
      arma::fill::zeros
  );

  arma::umat Parents(
      nVar,
      nSel,
      arma::fill::zeros
  );

  arma::uvec Best(
      nVar,
      arma::fill::zeros
  );

  arma::uvec Best_sim_phase1(
      nVar,
      arma::fill::zeros
  );


  arma::vec uProgeny(nPop, arma::fill::zeros), uParents(nSel, arma::fill::zeros);
  arma::vec simProgeny(nPop, arma::fill::zeros), simParents(nSel, arma::fill::zeros);
  arma::vec angleProgeny(nPop, arma::fill::zeros), angleParents(nSel, arma::fill::zeros);
  arma::vec lenProgeny(nPop, arma::fill::zeros), lenParents(nSel, arma::fill::zeros);
  arma::vec valProgeny(nPop, arma::fill::zeros), valParents(nSel, arma::fill::zeros);
  arma::uvec rankProgeny(nPop, arma::fill::zeros);

  double uBest_val=0.0, simBest=0.0, valBest=0.0, angleBest=0.0, lenBest=0.0;
  double uMax=0.0, uMin=0.0, simMax=0.0, simMin=0.0;
  double ufixedSum = arma::accu(ufixed);
  arma::uword currentRun = 0;

  // deterministic base seed (mix rd into 64-bit then fold)
  std::random_device rd;
  const uint64_t base_seed = (static_cast<uint64_t>(rd()) << 32) ^ static_cast<uint64_t>(rd());

  // contribution of fixed plan
  arma::vec xfixed =
    calcContr(
      fixedCrosses_mat,
      nInd,
      nCross
    );
  // ---------------------------------------------------------------------------
  // Initialize the additional constraints
  // ---------------------------------------------------------------------------
  
  // Transpose so each candidate cross is stored in one contiguous column.
  arma::mat scoreByConstraint =
    A.t();
  
  // Calculate additional constraint contributions from fixed crosses.
  arma::vec fixedExtra(
      A.n_cols,
      arma::fill::zeros
  );
  
  for (arma::uword f = 0; f < Afixed.n_rows; ++f) {
    
    for (arma::uword j = 0; j < A.n_cols; ++j) {
      
      fixedExtra(j) +=
        Afixed(f, j) / nCross;
    }
  }
  
  // Create one common constraint system.
  ConstraintSystem constraints{
    Crosses_mat,
    scoreByConstraint,
    minContr,
    maxContr,
    lower,
    upper,
    xfixed,
    fixedExtra,
    nCross,
    nInd,
    1e-12
  };
  // Max gain (no GA)
  arma::uvec uOrder =
    arma::sort_index(u, "descend");

  // Construct a high-gain reference plan satisfying
  // parental, QTL, and cross-type constraints.
  
  arma::uvec uBestIndex =
    build_constrained_umax_plan(
      u,
      uOrder,
      constraints,
      nVar,
      base_seed
    );

  {
    arma::vec x =
      calcContr(
        Crosses_mat.rows(uBestIndex),
        nInd,
        nCross
      );

    uMax =
      (
          arma::accu(u.elem(uBestIndex)) +
            ufixedSum
      ) /
        static_cast<double>(nCross);

    simMax =
      calc_sim(
        x + xfixed,
        G
      );
  }

  if (targetAngle < 1e-6) {
    
    arma::umat outCrossPlan =
      arma::join_cols(
        Crosses_mat.rows(uBestIndex),
        fixedCrosses_mat
      );
    
    ConstraintState bestState =
      constraints.initialize(uBestIndex);
    
    return Rcpp::List::create(
      Rcpp::Named("crossPlan") = outCrossPlan + 1,
      Rcpp::Named("uMax") = uMax,
      Rcpp::Named("simMax") = simMax,
      Rcpp::Named("parentContributions") = bestState.parents,
      Rcpp::Named("constraintValues") = bestState.extra
    );
  }

  // -------------------------
  // Phase 1: Optimize similarity (minimize sim)
  // -------------------------
  Rcpp::Rcout << "Optimize for Similarity" << std::endl << std::endl;

  ct_parallel_for(
    0,
    static_cast<int>(nPop),
    [&](int i)
    {
      const arma::uword ii =
        static_cast<arma::uword>(i);

      std::mt19937 gen(
          static_cast<uint32_t>(
            base_seed ^
              (
                  0x9E3779B97F4A7C15ULL +
                    static_cast<uint64_t>(i)
              )
          )
      );

      Progeny.col(ii) =
        sampleInt_std(
          nVar,
          potCross,
          gen
        );
    });

  repair_constraint_plans(
    Progeny,
    constraints,
    uBestIndex,
    base_seed ^ 0x1010101010101010ULL
  );

  ct_parallel_for(
    0,
    static_cast<int>(nPop),
    [&](int i)
    {
      const arma::uword ii =
        static_cast<arma::uword>(i);

      const arma::uvec idx =
        arma::conv_to<arma::uvec>::from(
          Progeny.col(ii)
        );

      const arma::vec x =
        calcContr(
          Crosses_mat.rows(idx),
          nInd,
          nCross
        );

      simProgeny(ii) =
        calc_sim(
          x + xfixed,
          G
        );
    });


  rankProgeny = arma::sort_index(simProgeny, "ascend");
  for (arma::uword i = 0; i < nSel; ++i) {
    Parents.col(i) = Progeny.col(rankProgeny(i));
    simParents(i)  = simProgeny(rankProgeny(i));
  }
  simBest =
    simParents(0);

  Best_sim_phase1 =
    arma::conv_to<arma::uvec>::from(
      Parents.col(0)
    );


  Rcpp::Rcout << "Gen  Similarity" << std::endl;

  for (arma::uword gen = 0; gen < maxGen; ++gen) {

  std::mt19937 rng_plan(
      static_cast<uint32_t>(
        base_seed ^
          (
              0xC3A5C85C97CB3127ULL +
                static_cast<uint64_t>(gen)
          )
      )
  );

    arma::umat crossPlan =
      sampHalfDialComb_std(
        nSel,
        nPop,
        rng_plan
      );
    ct_parallel_for(
      0,
      static_cast<int>(nPop),
      [&](int i)
      {
        const arma::uword ii =
          static_cast<arma::uword>(i);

        std::mt19937 gen_i(
            static_cast<uint32_t>(
              base_seed ^
                (
                    0xD2B74407B1CE6E93ULL +
                      static_cast<uint64_t>(gen) *
                      1315423911ULL +
                      static_cast<uint64_t>(i)
                )
            )
        );

        arma::uvec child =
          mate_uniform(
            arma::conv_to<arma::uvec>::from(
              Parents.col(
                crossPlan(0, ii)
              )
            ),
            arma::conv_to<arma::uvec>::from(
              Parents.col(
                crossPlan(1, ii)
              )
            ),
            gen_i,
            potCross
          );

        std::uniform_real_distribution<double>
          U01(0.0, 1.0);

        if (U01(gen_i) < probMut) {
          child =
            mutate(
              child,
              nMutate,
              potCross,
              gen_i
            );
        }

        Progeny.col(ii) = child;
      });
    repair_constraint_plans(
      Progeny,
      constraints,
      uBestIndex,
      base_seed ^
        (
            0x2020202020202020ULL +
              static_cast<uint64_t>(gen)
        )
    );
    
    ct_parallel_for(
      0,
      static_cast<int>(nPop),
      [&](int i)
      {
        const arma::uword ii =
          static_cast<arma::uword>(i);

        const arma::uvec idx =
          arma::conv_to<arma::uvec>::from(
            Progeny.col(ii)
          );

        const arma::vec x =
          calcContr(
            Crosses_mat.rows(idx),
            nInd,
            nCross
          );

        simProgeny(ii) =
          calc_sim(
            x + xfixed,
            G
          );
      });

    rankProgeny = arma::sort_index(simProgeny, "ascend");
    for (arma::uword i = 0; i < nSel; ++i) {
      Parents.col(i) = Progeny.col(rankProgeny(i));
      simParents(i)  = simProgeny(rankProgeny(i));
    }

    if (simParents(0) < simBest) {
      simBest =
        simParents(0);

      Best_sim_phase1 =
        arma::conv_to<arma::uvec>::from(
          Parents.col(0)
        );

      currentRun = 0;

    } else {
      ++currentRun;
    }

    if (gen % 10 == 0) Rcpp::Rcout << gen << "  " << simBest << std::endl;
    if (currentRun >= maxRun) break;
  }

  simMin =
    simBest;

  uMin =
    (
        arma::accu(
          u.elem(Best_sim_phase1)
        ) +
          ufixedSum
    ) /
      static_cast<double>(nCross);


  // -------------------------
  // Phase 2: Optimize crossing plan (angle/length w.r.t. targetAngle)
  // -------------------------
  Rcpp::Rcout << std::endl << std::endl << "Optimize for Crossing Plan" << std::endl << std::endl;

  ct_parallel_for(
    0,
    static_cast<int>(nPop),
    [&](int i)
    {
      const arma::uword ii =
        static_cast<arma::uword>(i);

      std::mt19937 gen_i(
          static_cast<uint32_t>(
            base_seed ^
              (
                  0xA24BAED4963EE407ULL +
                    static_cast<uint64_t>(i)
              )
          )
      );

      Progeny.col(ii) =
        sampleInt_std(
          nVar,
          potCross,
          gen_i
        );
    });

  if (nPop > 0) {
    Progeny.col(0) = uBestIndex;
  }

  if (nPop > 1) {
    Progeny.col(1) = Best_sim_phase1;
  }

  repair_constraint_plans(
    Progeny,
    constraints,
    uBestIndex,
    base_seed ^ 0x3030303030303030ULL
  );

  ct_parallel_for(
    0,
    static_cast<int>(nPop),
    [&](int i)
    {
      const arma::uword ii =
        static_cast<arma::uword>(i);

      const arma::uvec idx =
        arma::conv_to<arma::uvec>::from(
          Progeny.col(ii)
        );

      const arma::vec x =
        calcContr(
          Crosses_mat.rows(idx),
          nInd,
          nCross
        );

      simProgeny(ii) =
        calc_sim(
          x + xfixed,
          G
        );

      uProgeny(ii) =
        (
            arma::accu(
              u.elem(idx)
            ) +
              ufixedSum
        ) /
          static_cast<double>(nCross);

      calcVec(
        angleProgeny(ii),
        lenProgeny(ii),
        uProgeny(ii),
        simProgeny(ii),
        uMax,
        simMax,
        uMin,
        simMin
      );

      valProgeny(ii) =
        lenProgeny(ii) -
        anglePenalty *
        std::abs(
          angleProgeny(ii) -
            targetAngle
        );
    });



  rankProgeny = arma::sort_index(valProgeny, "descend");
  for (arma::uword i = 0; i < nSel; ++i) {
    Parents.col(i) = Progeny.col(rankProgeny(i));
    uParents(i)    = uProgeny(rankProgeny(i));
    simParents(i)  = simProgeny(rankProgeny(i));
    angleParents(i)= angleProgeny(rankProgeny(i));
    lenParents(i)  = lenProgeny(rankProgeny(i));
    valParents(i)  = valProgeny(rankProgeny(i));
  }

  valBest   = valParents(0);
  angleBest = angleParents(0);
  lenBest   = lenParents(0);
  uBest_val = uParents(0);
  simBest   = simParents(0);
  Best      = Parents.col(0);

  Rcpp::Rcout << "Gen  acc_cross_value  Similarity  Angle  Length  Value" << std::endl;
  currentRun = 0;

  for (arma::uword gen = 0; gen < maxGen; ++gen) {

    std::mt19937 rng_plan(
        static_cast<uint32_t>(
          base_seed ^
            (
                0x9E3779B185EBCA87ULL +
                  static_cast<uint64_t>(gen)
            )
        )
    );

    arma::umat crossPlan =
      sampHalfDialComb_std(
        nSel,
        nPop,
        rng_plan
      );

    ct_parallel_for(
      0,
      static_cast<int>(nPop),
      [&](int i)
      {
        const arma::uword ii =
          static_cast<arma::uword>(i);

        std::mt19937 gen_i(
            static_cast<uint32_t>(
              base_seed ^
                (
                    0xC949D7C7509E6557ULL +
                      static_cast<uint64_t>(gen) *
                      40503ULL +
                      static_cast<uint64_t>(i)
                )
            )
        );

        arma::uvec child =
          mate_uniform(
            arma::conv_to<arma::uvec>::from(
              Parents.col(
                crossPlan(0, ii)
              )
            ),
            arma::conv_to<arma::uvec>::from(
              Parents.col(
                crossPlan(1, ii)
              )
            ),
            gen_i,
            potCross
          );

        std::uniform_real_distribution<double>
          U01(0.0, 1.0);

        if (U01(gen_i) < probMut) {
          child =
            mutate(
              child,
              nMutate,
              potCross,
              gen_i
            );
        }

        Progeny.col(ii) = child;
      });

    repair_constraint_plans(
      Progeny,
      constraints,
      uBestIndex,
      base_seed ^
        (
            0x4040404040404040ULL +
              static_cast<uint64_t>(gen)
        )
    );

    ct_parallel_for(
      0,
      static_cast<int>(nPop),
      [&](int i)
      {
        const arma::uword ii =
          static_cast<arma::uword>(i);

        const arma::uvec idx =
          arma::conv_to<arma::uvec>::from(
            Progeny.col(ii)
          );

        const arma::vec x =
          calcContr(
            Crosses_mat.rows(idx),
            nInd,
            nCross
          );

        simProgeny(ii) =
          calc_sim(
            x + xfixed,
            G
          );

        uProgeny(ii) =
          (
              arma::accu(
                u.elem(idx)
              ) +
                ufixedSum
          ) /
            static_cast<double>(nCross);

        calcVec(
          angleProgeny(ii),
          lenProgeny(ii),
          uProgeny(ii),
          simProgeny(ii),
          uMax,
          simMax,
          uMin,
          simMin
        );

        valProgeny(ii) =
          lenProgeny(ii) -
          anglePenalty *
          std::abs(
            angleProgeny(ii) -
              targetAngle
          );
      });



    rankProgeny = arma::sort_index(valProgeny, "descend");
    for (arma::uword i = 0; i < nSel; ++i) {
      Parents.col(i)   = Progeny.col(rankProgeny(i));
      uParents(i)      = uProgeny(rankProgeny(i));
      simParents(i)    = simProgeny(rankProgeny(i));
      angleParents(i)  = angleProgeny(rankProgeny(i));
      lenParents(i)    = lenProgeny(rankProgeny(i));
      valParents(i)    = valProgeny(rankProgeny(i));
    }

    if (valParents(0) > valBest) {
      valBest   = valParents(0);
      angleBest = angleParents(0);
      lenBest   = lenParents(0);
      uBest_val = uParents(0);
      simBest   = simParents(0);
      Best      = Parents.col(0);
      currentRun = 0;
    } else {
      ++currentRun;
    }

    if (gen % 10 == 0) {
      Rcpp::Rcout << gen << "  " << uBest_val << "  " << simBest << "  "
                  << angleBest << "  " << lenBest << "  " << valBest << std::endl;
    }
    if (currentRun >= maxRun) break;
  }

  // ---------------------------------------------------------------------------
  // Verify all constraint classes in the final mating plan
  // ---------------------------------------------------------------------------
  
  ConstraintState finalState =
    constraints.initialize(Best);
  
  if (finalState.violation > 1e-10) {
    
    Rcpp::stop(
      "Internal error: final mating plan violates constraints."
    );
  }
  
  arma::umat outCrossPlan =
    arma::join_cols(
      Crosses_mat.rows(Best),
      fixedCrosses_mat
    );


  return Rcpp::List::create(
    Rcpp::Named("crossPlan") = outCrossPlan + 1,
    Rcpp::Named("uMax") = uMax,
    Rcpp::Named("uMin") = uMin,
    Rcpp::Named("simMax") = simMax,
    Rcpp::Named("simMin") = simMin,
    Rcpp::Named("uBest") = uBest_val,
    Rcpp::Named("simBest") = simBest,
    Rcpp::Named("angleBest") = angleBest,
    Rcpp::Named("lenBest") = lenBest,
    
    // Achieved parental contributions.
    Rcpp::Named("parentContributions") =
      finalState.parents,
      
      // Achieved QTL frequencies and cross-type proportions.
    Rcpp::Named("constraintValues") =
      finalState.extra
  );
}
