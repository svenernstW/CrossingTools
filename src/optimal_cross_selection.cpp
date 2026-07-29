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
static arma::uvec sampleInt_std(arma::uword n, arma::uword N, std::mt19937& gen) {
  if (n > N) n = N;
  std::vector<arma::uword> pool(N);
  std::iota(pool.begin(), pool.end(), 0);
  // partial Fisher–Yates: bring n items to front
  for (arma::uword i = 0; i < n; ++i) {
    std::uniform_int_distribution<arma::uword> dis(i, N - 1);
    arma::uword j = dis(gen);
    std::swap(pool[i], pool[j]);
  }
  std::vector<arma::uword> out(pool.begin(), pool.begin() + n);
  std::sort(out.begin(), out.end());
  return arma::uvec(out);
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

static inline arma::uvec build_constrained_umax_plan(
    const arma::vec& u,
    const arma::uvec& uOrder,       // cross indices sorted by descending u
    const arma::umat& Crosses_mat,  // potCross x 2
    arma::uword nVar,
    arma::uword nCross,
    arma::uword nInd,
    const arma::vec& xfixed,
    const arma::vec& minContr,
    const arma::vec& maxContr,
    double tol = 1e-12
) {
  const arma::uword potCross = Crosses_mat.n_rows;
  
  if (nCross == 0) {
    Rcpp::stop("nCross must be greater than zero.");
  }
  
  if (nVar > potCross) {
    Rcpp::stop(
      "Cannot construct uMax plan: nVar exceeds the number of candidate crosses."
    );
  }
  
  if (minContr.n_elem != nInd || maxContr.n_elem != nInd) {
    Rcpp::stop(
      "minContr and maxContr must contain one value per parent."
    );
  }
  
  auto value_violation = [&](
    arma::uword parent,
    double value
  ) -> double {
    double out = 0.0;
    
    if (value < minContr(parent) - tol) {
      out += minContr(parent) - value;
    }
    
    if (value > maxContr(parent) + tol) {
      out += value - maxContr(parent);
    }
    
    return out;
  };
  
  // Handle a plan containing only fixed crosses.
  if (nVar == 0) {
    double violation = 0.0;
    
    for (arma::uword p = 0; p < nInd; ++p) {
      violation += value_violation(p, xfixed(p));
    }
    
    if (violation > tol) {
      Rcpp::stop(
        "The contribution constraints cannot be satisfied by the fixed crosses."
      );
    }
    
    return arma::uvec();
  }
  
  // Start with the unconstrained top-u plan.
  arma::uvec plan = uOrder.head(nVar);
  
  std::vector<unsigned char> selected(
      static_cast<size_t>(potCross),
      0
  );
  
  for (arma::uword k = 0; k < nVar; ++k) {
    selected[static_cast<size_t>(plan(k))] = 1;
  }
  
  arma::vec contr =
    calcContr(
      Crosses_mat.rows(plan),
      nInd,
      nCross
    ) + xfixed;
  
  const double step =
    1.0 / (2.0 * static_cast<double>(nCross));
  
  auto total_violation = [&]() -> double {
    double out = 0.0;
    
    for (arma::uword p = 0; p < nInd; ++p) {
      out += value_violation(p, contr(p));
    }
    
    return out;
  };
  
  /*
   * Compute constraint violation after replacing plan(pos)
   * with newCross. The contribution vector is restored before
   * returning.
   */
  auto violation_after_swap = [&](
    arma::uword pos,
    arma::uword newCross,
    double currentViolation
  ) -> double {
    const arma::uword oldCross = plan(pos);
    
    const arma::uword oldA = Crosses_mat(oldCross, 0);
    const arma::uword oldB = Crosses_mat(oldCross, 1);
    const arma::uword newA = Crosses_mat(newCross, 0);
    const arma::uword newB = Crosses_mat(newCross, 1);
    
    arma::uword affected[4];
    int nAffected = 0;
    
    auto add_affected = [&](arma::uword parent) {
      for (int j = 0; j < nAffected; ++j) {
        if (affected[j] == parent) return;
      }
      
      affected[nAffected++] = parent;
    };
    
    add_affected(oldA);
    add_affected(oldB);
    add_affected(newA);
    add_affected(newB);
    
    double before = 0.0;
    
    for (int j = 0; j < nAffected; ++j) {
      const arma::uword p = affected[j];
      before += value_violation(p, contr(p));
    }
    
    contr(oldA) -= step;
    contr(oldB) -= step;
    contr(newA) += step;
    contr(newB) += step;
    
    double after = 0.0;
    
    for (int j = 0; j < nAffected; ++j) {
      const arma::uword p = affected[j];
      after += value_violation(p, contr(p));
    }
    
    // Undo tentative change.
    contr(oldA) += step;
    contr(oldB) += step;
    contr(newA) -= step;
    contr(newB) -= step;
    
    return std::max(
      0.0,
      currentViolation - before + after
    );
  };
  
  auto apply_swap = [&](
    arma::uword pos,
    arma::uword newCross
  ) {
    const arma::uword oldCross = plan(pos);
    
    const arma::uword oldA = Crosses_mat(oldCross, 0);
    const arma::uword oldB = Crosses_mat(oldCross, 1);
    const arma::uword newA = Crosses_mat(newCross, 0);
    const arma::uword newB = Crosses_mat(newCross, 1);
    
    contr(oldA) -= step;
    contr(oldB) -= step;
    contr(newA) += step;
    contr(newB) += step;
    
    selected[static_cast<size_t>(oldCross)] = 0;
    selected[static_cast<size_t>(newCross)] = 1;
    
    plan(pos) = newCross;
  };
  
  /*
   * Phase 1:
   * Repair the top-u plan.
   */
  double currentViolation = total_violation();
  
  const arma::uword maxRepairIterations =
    std::max<arma::uword>(100, 10 * nVar);
  
  for (
      arma::uword iter = 0;
      iter < maxRepairIterations && currentViolation > tol;
      ++iter
  ) {
    std::vector<unsigned char> under(
        static_cast<size_t>(nInd),
        0
    );
    
    bool anyUnder = false;
    
    for (arma::uword p = 0; p < nInd; ++p) {
      if (contr(p) < minContr(p) - tol) {
        under[static_cast<size_t>(p)] = 1;
        anyUnder = true;
      }
    }
    
    // Prefer removing the lowest-u selected cross.
    std::vector<arma::uword> removeOrder(
        static_cast<size_t>(nVar)
    );
    
    std::iota(
      removeOrder.begin(),
      removeOrder.end(),
      static_cast<arma::uword>(0)
    );
    
    std::sort(
      removeOrder.begin(),
      removeOrder.end(),
      [&](arma::uword a, arma::uword b) {
        return u(plan(a)) < u(plan(b));
      }
    );
    
    bool foundSwap = false;
    arma::uword selectedPos = 0;
    arma::uword selectedNewCross = 0;
    double selectedViolation = currentViolation;
    
    /*
     * New crosses are tested from highest to lowest u.
     * For each new cross, old crosses are tested from
     * lowest to highest u.
     */
    for (
        arma::uword ord = 0;
        ord < uOrder.n_elem && !foundSwap;
        ++ord
    ) {
      const arma::uword newCross = uOrder(ord);
      
      if (selected[static_cast<size_t>(newCross)]) {
        continue;
      }
      
      const arma::uword newA = Crosses_mat(newCross, 0);
      const arma::uword newB = Crosses_mat(newCross, 1);
      
      /*
       * When a minimum is violated, require the new cross
       * to increase at least one currently deficient parent.
       */
      if (
          anyUnder &&
            !under[static_cast<size_t>(newA)] &&
            !under[static_cast<size_t>(newB)]
      ) {
        continue;
      }
      
      for (arma::uword pos : removeOrder) {
        const double newViolation =
          violation_after_swap(
            pos,
            newCross,
            currentViolation
          );
        
        if (newViolation < currentViolation - tol) {
          selectedPos = pos;
          selectedNewCross = newCross;
          selectedViolation = newViolation;
          foundSwap = true;
          break;
        }
      }
    }
    
    if (!foundSwap) {
      Rcpp::stop(
        "Could not construct a contribution-feasible uMax plan "
        "using single-cross replacements. The constraints may be "
        "infeasible or may require a multi-cross replacement."
      );
    }
    
    apply_swap(
      selectedPos,
      selectedNewCross
    );
    
    selectedViolation = total_violation();
    currentViolation = selectedViolation;
  }
  
  currentViolation = total_violation();
  
  if (currentViolation > tol) {
    Rcpp::stop(
      "Contribution repair for the uMax plan did not converge."
    );
  }
  
  /*
   * Phase 2:
   * Feasible 1-swap local search.
   *
   * Replace a selected cross with a higher-u unselected cross
   * whenever all contribution constraints remain satisfied.
   */
  bool improved = true;
  
  while (improved) {
    improved = false;
    
    std::vector<arma::uword> removeOrder(
        static_cast<size_t>(nVar)
    );
    
    std::iota(
      removeOrder.begin(),
      removeOrder.end(),
      static_cast<arma::uword>(0)
    );
    
    std::sort(
      removeOrder.begin(),
      removeOrder.end(),
      [&](arma::uword a, arma::uword b) {
        return u(plan(a)) < u(plan(b));
      }
    );
    
    const double lowestSelectedU =
      u(plan(removeOrder.front()));
    
    for (
        arma::uword ord = 0;
        ord < uOrder.n_elem && !improved;
        ++ord
    ) {
      const arma::uword newCross = uOrder(ord);
      
      if (selected[static_cast<size_t>(newCross)]) {
        continue;
      }
      
      /*
       * Since uOrder is descending, once the new cross is
       * no better than the lowest-u selected cross, no later
       * candidate can improve the objective.
       */
      if (u(newCross) <= lowestSelectedU + tol) {
        break;
      }
      
      for (arma::uword pos : removeOrder) {
        const arma::uword oldCross = plan(pos);
        
        if (u(newCross) <= u(oldCross) + tol) {
          break;
        }
        
        const double newViolation =
          violation_after_swap(
            pos,
            newCross,
            0.0
          );
        
        if (newViolation <= tol) {
          apply_swap(pos, newCross);
          improved = true;
          break;
        }
      }
    }
  }
  
  if (total_violation() > tol) {
    Rcpp::stop(
      "Internal error: constrained uMax plan is not feasible."
    );
  }
  
  return plan;
}

static inline void fix_contribution_plans(
    arma::umat& Plans,              // nVar x nPlans
    const arma::umat& Crosses_mat,  // potCross x 2, 0-based parents
    const arma::vec& minContr,      // fractions, e.g. 0.07
    const arma::vec& maxContr,      // fractions or Inf
    const arma::vec& xfixed,        // already normalized by nCross
    arma::uword nCross,
    arma::uword nInd,
    uint64_t base_seed,
    int maxIter = 5000,
    int maxCandidateTries = 1000000,
    double tol = 1e-12
) {
  const arma::uword nVar     = Plans.n_rows;
  const arma::uword nPlans   = Plans.n_cols;
  const arma::uword potCross = Crosses_mat.n_rows;
  
  if (nPlans == 0) return;
  
  // Completely unconstrained case:
  // min = 0 and max = Inf for every individual.
  bool noConstraints = true;
  
  for (arma::uword p = 0; p < nInd; ++p) {
    if (minContr(p) > tol || std::isfinite(maxContr(p))) {
      noConstraints = false;
      break;
    }
  }
  
  if (noConstraints) return;
  
  auto value_violation = [&](arma::uword p, double value) -> double {
    double v = 0.0;
    
    if (value < minContr(p) - tol) {
      v += minContr(p) - value;
    }
    
    if (value > maxContr(p) + tol) {
      v += value - maxContr(p);
    }
    
    return v;
  };
  
  // No variable crosses: only verify the fixed plan.
  if (nVar == 0) {
    double v = 0.0;
    
    for (arma::uword p = 0; p < nInd; ++p) {
      v += value_violation(p, xfixed(p));
    }
    
    if (v > tol) {
      Rcpp::stop(
        "Contribution constraints cannot be satisfied because "
        "the plan contains only fixed crosses."
      );
    }
    
    return;
  }
  
  if (potCross == 0) {
    Rcpp::stop("No candidate crosses are available.");
  }
  
  const double step =
    1.0 / (2.0 * static_cast<double>(nCross));
  
  // Written only once per plan, so separate int elements are thread-safe.
  std::vector<int> failed(
      static_cast<size_t>(nPlans), 0
  );
  
  ct_parallel_for(0, static_cast<int>(nPlans), [&](int ii) {
    const arma::uword planIndex =
      static_cast<arma::uword>(ii);
    
    std::mt19937 gen(
        static_cast<uint32_t>(
          base_seed ^
            (0xA13B77C19E3779B9ULL +
            static_cast<uint64_t>(planIndex))
        )
    );
    
    arma::uvec plan =
      arma::conv_to<arma::uvec>::from(
        Plans.col(planIndex)
      );
    
    // Only stores the selected crosses, rather than allocating
    // potCross bytes for every plan.
    std::unordered_set<arma::uword> inPlan;
    inPlan.reserve(
      static_cast<size_t>(nVar) * 2 + 1
    );
    
    for (arma::uword k = 0; k < nVar; ++k) {
      inPlan.insert(plan(k));
    }
    
    // Both variable and fixed contributions use nCross as denominator.
    arma::vec contr =
      calcContr(
        Crosses_mat.rows(plan),
        nInd,
        nCross
      ) + xfixed;
    
    auto total_violation = [&]() -> double {
      double v = 0.0;
      
      for (arma::uword p = 0; p < nInd; ++p) {
        v += value_violation(p, contr(p));
      }
      
      return v;
    };
    
    double currentV = total_violation();
    
    for (
        int iter = 0;
        iter < maxIter && currentV > tol;
        ++iter
    ) {
      std::vector<unsigned char> underFlag(nInd, 0);
      std::vector<unsigned char> overFlag(nInd, 0);
      
      bool anyUnder = false;
      bool anyOver  = false;
      
      for (arma::uword p = 0; p < nInd; ++p) {
        if (contr(p) < minContr(p) - tol) {
          underFlag[p] = 1;
          anyUnder = true;
        }
        
        if (contr(p) > maxContr(p) + tol) {
          overFlag[p] = 1;
          anyOver = true;
        }
      }
      
      if (!anyUnder && !anyOver) {
        currentV = 0.0;
        break;
      }
      
      std::vector<arma::uword> removePos;
      removePos.reserve(static_cast<size_t>(nVar));
      
      if (anyUnder) {
        /*
         * Best case:
         * remove a cross containing an overrepresented parent
         * but no underrepresented parent.
         */
        if (anyOver) {
          for (arma::uword k = 0; k < nVar; ++k) {
            const arma::uword cr = plan(k);
            const arma::uword a  = Crosses_mat(cr, 0);
            const arma::uword b  = Crosses_mat(cr, 1);
            
            const bool containsUnder =
              underFlag[a] || underFlag[b];
            
            const bool containsOver =
              overFlag[a] || overFlag[b];
            
            if (!containsUnder && containsOver) {
              removePos.push_back(k);
            }
          }
        }
        
        /*
         * Second preference:
         * remove any cross that does not contain an
         * underrepresented parent.
         */
        if (removePos.empty()) {
          for (arma::uword k = 0; k < nVar; ++k) {
            const arma::uword cr = plan(k);
            const arma::uword a  = Crosses_mat(cr, 0);
            const arma::uword b  = Crosses_mat(cr, 1);
            
            if (!underFlag[a] && !underFlag[b]) {
              removePos.push_back(k);
            }
          }
        }
        
        /*
         * Fallback:
         * allow any removal. The swap is accepted only
         * if total violation decreases.
         */
        if (removePos.empty()) {
          for (arma::uword k = 0; k < nVar; ++k) {
            removePos.push_back(k);
          }
        }
        
      } else {
        // Only maximum constraints remain.
        for (arma::uword k = 0; k < nVar; ++k) {
          const arma::uword cr = plan(k);
          const arma::uword a  = Crosses_mat(cr, 0);
          const arma::uword b  = Crosses_mat(cr, 1);
          
          if (overFlag[a] || overFlag[b]) {
            removePos.push_back(k);
          }
        }
      }
      
      if (removePos.empty()) break;
      
      std::uniform_int_distribution<size_t> disRemove(
          0,
          removePos.size() - 1
      );
      
      std::uniform_int_distribution<arma::uword> disCross(
          0,
          potCross - 1
      );
      
      bool accepted = false;
      
      for (
          int tr = 0;
          tr < maxCandidateTries && !accepted;
          ++tr
      ) {
        const arma::uword pos =
          removePos[disRemove(gen)];
        
        const arma::uword oldCr = plan(pos);
        const arma::uword newCr = disCross(gen);
        
        // Crosses within one plan must stay unique.
        if (inPlan.find(newCr) != inPlan.end()) {
          continue;
        }
        
        const arma::uword oldA =
          Crosses_mat(oldCr, 0);
        
        const arma::uword oldB =
          Crosses_mat(oldCr, 1);
        
        const arma::uword newA =
          Crosses_mat(newCr, 0);
        
        const arma::uword newB =
          Crosses_mat(newCr, 1);
        
        // When minimum constraints are violated, the added
        // cross must contain at least one underrepresented parent.
        if (
            anyUnder &&
              !underFlag[newA] &&
              !underFlag[newB]
        ) {
          continue;
        }
        
        // Collect the unique affected parents.
        arma::uword affected[4];
        int nAffected = 0;
        
        auto add_affected = [&](arma::uword p) {
          for (int j = 0; j < nAffected; ++j) {
            if (affected[j] == p) return;
          }
          
          affected[nAffected++] = p;
        };
        
        add_affected(oldA);
        add_affected(oldB);
        add_affected(newA);
        add_affected(newB);
        
        double affectedBefore = 0.0;
        
        for (int j = 0; j < nAffected; ++j) {
          const arma::uword p = affected[j];
          
          affectedBefore +=
            value_violation(p, contr(p));
        }
        
        // Apply the tentative swap.
        contr(oldA) -= step;
        contr(oldB) -= step;
        contr(newA) += step;
        contr(newB) += step;
        
        double affectedAfter = 0.0;
        
        for (int j = 0; j < nAffected; ++j) {
          const arma::uword p = affected[j];
          
          affectedAfter +=
            value_violation(p, contr(p));
        }
        
        double newV =
          currentV -
          affectedBefore +
          affectedAfter;
        
        newV = std::max(0.0, newV);
        
        if (newV < currentV - tol) {
          inPlan.erase(oldCr);
          inPlan.insert(newCr);
          
          plan(pos) = newCr;
          currentV = newV;
          accepted = true;
          
        } else {
          // Undo rejected swap.
          contr(oldA) += step;
          contr(oldB) += step;
          contr(newA) -= step;
          contr(newB) -= step;
        }
      }
      
      if (!accepted) break;
    }
    
    // Full recalculation to avoid accumulated floating-point drift.
    currentV = total_violation();
    
    if (currentV > tol) {
      failed[static_cast<size_t>(planIndex)] = 1;
    }
    
    Plans.col(planIndex) = plan;
  });
  
  size_t nFailed = 0;
  size_t firstFailed = 0;
  
  for (size_t i = 0; i < failed.size(); ++i) {
    if (failed[i]) {
      if (nFailed == 0) firstFailed = i;
      ++nFailed;
    }
  }
  
  if (nFailed > 0) {
    std::ostringstream msg;
    
    msg
    << "Contribution repair failed for "
    << nFailed
    << " plan(s). First failed plan: "
    << (firstFailed + 1)
    << ". The constraints may be infeasible or too restrictive.";
    
    Rcpp::stop(msg.str());
  }
}

// Normalize to unit square, compute angle and length toward (gain high, similarity low)
static void calcVec(double& angle, double& length, double u, double sim,
             const double& uMax, const double& simMax,
             const double& uMin, const double& simMin) {
  u   = (u   - uMin) / (uMax - uMin);
  sim = (simMax - sim) / (simMax - simMin);
  length = std::sqrt(u * u + sim * sim);
  angle  = std::acos(u / std::max(length, 1e-16));
  if (u < 0) length = -length;
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
  
  if (n == 0 || potCross == 0) {
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
    double probMut=0.01,
    arma::uword nMutate=2,
    arma::uword nSel=500,
    arma::uword nPop=10000,
    arma::uword maxGen=1000,
    arma::uword maxRun=100,
    double anglePenalty=0.5,
    int nThreads=4){
  // portable threading
  ct_set_threads(nThreads);

  const arma::uword potCross  = Crosses.nrow();
  const arma::uword fixedRows = fixedCrosses.nrow();
  const arma::uword nVar      = nCross - fixedRows;

  arma::umat Crosses_mat      = as<arma::umat>(Crosses);
  arma::umat fixedCrosses_mat = as<arma::umat>(fixedCrosses);
  Crosses_mat      -= 1; // to 0-based
  fixedCrosses_mat -= 1;

  arma::uword nInd = G.n_cols;

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
  
  // Max gain (no GA)
  arma::uvec uOrder =
    arma::sort_index(u, "descend");
  
  bool hasContributionConstraints = false;
  
  for (arma::uword p = 0; p < nInd; ++p) {
    const bool bindingMinimum =
      minContr(p) > 1e-12;
    
    const bool bindingMaximum =
      std::isfinite(maxContr(p)) &&
      maxContr(p) < 1.0 - 1e-12;
    
    if (bindingMinimum || bindingMaximum) {
      hasContributionConstraints = true;
      break;
    }
  }
  
  arma::uvec uBestIndex;
  
  if (hasContributionConstraints) {
    uBestIndex =
      build_constrained_umax_plan(
        u,
        uOrder,
        Crosses_mat,
        nVar,
        nCross,
        nInd,
        xfixed,
        minContr,
        maxContr
      );
  } else {
    uBestIndex = uOrder.head(nVar);
  }
  
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
    arma::umat outCrossPlan = arma::join_cols(Crosses_mat.rows(uBestIndex), fixedCrosses_mat);
    return Rcpp::List::create(
      Rcpp::Named("crossPlan") = outCrossPlan + 1,
      Rcpp::Named("uMax")      = uMax,
      Rcpp::Named("simMax")    = simMax
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
  
  fix_contribution_plans(
    Progeny,
    Crosses_mat,
    minContr,
    maxContr,
    xfixed,
    nCross,
    nInd,
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
    fix_contribution_plans(
      Progeny,
      Crosses_mat,
      minContr,
      maxContr,
      xfixed,
      nCross,
      nInd,
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
  
  fix_contribution_plans(
    Progeny,
    Crosses_mat,
    minContr,
    maxContr,
    xfixed,
    nCross,
    nInd,
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
    
    fix_contribution_plans(
      Progeny,
      Crosses_mat,
      minContr,
      maxContr,
      xfixed,
      nCross,
      nInd,
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

  arma::vec finalContribution =
    calcContr(
      Crosses_mat.rows(Best),
      nInd,
      nCross
    ) + xfixed;
  
  double maxContributionViolation =
    0.0;
  
  for (arma::uword p = 0; p < nInd; ++p) {
    if (
        finalContribution(p) <
          minContr(p)
    ) {
      maxContributionViolation =
        std::max(
          maxContributionViolation,
          minContr(p) -
            finalContribution(p)
        );
    }
    
    if (
        finalContribution(p) >
      maxContr(p)
    ) {
      maxContributionViolation =
        std::max(
          maxContributionViolation,
          finalContribution(p) -
            maxContr(p)
        );
    }
  }
  
  if (maxContributionViolation > 1e-10) {
    Rcpp::stop(
      "Internal error: final crossing plan violates "
      "the contribution constraints."
    );
  }
  
  Best =
    arma::sort(Best);
  
  arma::umat outCrossPlan =
    arma::join_cols(
      Crosses_mat.rows(Best),
      fixedCrosses_mat
    );
  
  
  return Rcpp::List::create(
    Rcpp::Named("crossPlan") = outCrossPlan + 1,
    Rcpp::Named("uMax")      = uMax,
    Rcpp::Named("uMin")      = uMin,
    Rcpp::Named("simMax")    = simMax,
    Rcpp::Named("simMin")    = simMin,
    Rcpp::Named("uBest")     = uBest_val,
    Rcpp::Named("simBest")   = simBest,
    Rcpp::Named("angleBest") = angleBest,
    Rcpp::Named("lenBest")   = lenBest
  );
}
