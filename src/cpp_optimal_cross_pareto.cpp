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

#include "parallel_backend.h"  // OpenMP when available, else RcppParallel

using namespace Rcpp;

// ------------------------------
// Helpers
// ------------------------------
// Parental contributions for a crossing plan (indices are 0-based)
static arma::vec calcContr(const arma::umat& crosses, const arma::uword& nInd) {
  const arma::uword M = crosses.n_rows;
  double val = 1.0 / (2.0 * static_cast<double>(M));
  arma::vec x(nInd, arma::fill::zeros);
  for (arma::uword r = 0; r < M; ++r) {
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
static inline arma::uvec sampleInt_std(arma::uword n, arma::uword N, std::mt19937& gen) {
  if (n > N) n = N;
  std::vector<arma::uword> pool(N);
  std::iota(pool.begin(), pool.end(), 0);

  for (arma::uword i = 0; i < n; ++i) {
    std::uniform_int_distribution<arma::uword> dis(i, N - 1);
    arma::uword j = dis(gen);
    std::swap(pool[i], pool[j]);
  }

  std::vector<arma::uword> out(pool.begin(), pool.begin() + n);
  std::sort(out.begin(), out.end());
  return arma::uvec(out);
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

// uniform set-based crossover
static inline arma::uvec mate_uniform(const arma::uvec& parent1, const arma::uvec& parent2,
                                      std::mt19937& gen, arma::uword potCross)
{
  arma::uword n = parent1.n_elem;
  arma::uvec offspring(n, arma::fill::zeros);
  std::vector<char> present(potCross, 0);
  offspring.fill(potCross);

  std::uniform_real_distribution<> dis(0.0, 1.0);
  for (arma::uword i = 0; i < n; ++i) {
    bool take_from_p1 = dis(gen) < 0.5;
    arma::uword gene = take_from_p1 ? parent1(i) : parent2(i);
    if (!present[gene]) {
      offspring(i) = gene; present[gene] = 1;
    } else {
      gene = take_from_p1 ? parent2(i) : parent1(i);
      if (!present[gene]) {
        offspring(i) = gene; present[gene] = 1;
      } else {
        offspring(i) = potCross;
      }
    }
  }

  for (arma::uword i = 0; i < n; ++i) if (offspring(i) == potCross) {
    arma::uword g1 = parent1(i);
    if (!present[g1]) { offspring(i) = g1; present[g1] = 1; continue; }
    arma::uword g2 = parent2(i);
    if (!present[g2]) { offspring(i) = g2; present[g2] = 1; continue; }
    for (arma::uword g = 0; g < potCross; ++g) if (!present[g]) {
      offspring(i) = g; present[g] = 1; break;
    }
  }

  for (arma::uword i = 0; i < n; ++i) if (offspring(i) >= potCross) {
    for (arma::uword g = 0; g < potCross; ++g) if (!present[g]) {
      offspring(i) = g; present[g] = 1; break;
    }
  }

  return offspring;
}

// mutation: replace nMutate positions with genes not in plan
static inline arma::uvec mutate(const arma::uvec& crosses, arma::uword nMutate,
                                arma::uword potCross, std::mt19937& gen)
{
  arma::uword n = crosses.n_elem;
  if (n == 0 || potCross == 0) return crosses;

  arma::uword k = std::min<arma::uword>(nMutate, n);
  if (k == 0) return crosses;

  arma::uvec out = crosses;

  std::vector<char> present(potCross, 0);
  for (arma::uword i = 0; i < n; ++i) {
    arma::uword g = out(i);
    if (g < potCross) present[g] = 1;
  }

  arma::uvec pos = sampleInt_std(k, n, gen);
  std::uniform_int_distribution<arma::uword> dis(0, potCross - 1);

  for (arma::uword t = 0; t < pos.n_elem; ++t) {
    arma::uword p = pos(t);
    arma::uword oldg = out(p);
    if (oldg < potCross) present[oldg] = 0;

    arma::uword newg = potCross;
    for (int tries = 0; tries < 64; ++tries) {
      arma::uword cand = dis(gen);
      if (!present[cand]) { newg = cand; break; }
    }
    if (newg == potCross) {
      for (arma::uword cand = 0; cand < potCross; ++cand) if (!present[cand]) { newg = cand; break; }
    }

    out(p) = newg;
    present[newg] = 1;
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
static inline void enforce_unique_population(
    arma::umat& Progeny,
    arma::uword nVar,
    arma::uword potCross,
    uint64_t base_seed)
{
  std::unordered_set<std::string> seen;
  seen.reserve((size_t)Progeny.n_cols * 2);

  for (arma::uword i = 0; i < Progeny.n_cols; ++i) {
    int tries = 0;
    for (;;) {
      arma::uvec v  = arma::conv_to<arma::uvec>::from(Progeny.col(i));
      if (v.n_elem != nVar) {
        // hard reset if corrupted
        std::mt19937 gen((uint32_t)(base_seed ^ (0xABCDEF0123456789ULL + (uint64_t)i)));
        Progeny.col(i) = sampleInt_std(nVar, potCross, gen);
        v = arma::conv_to<arma::uvec>::from(Progeny.col(i));
      }

      arma::uvec vs = arma::sort(v);
      std::string key = plan_key_sorted(vs);
      if (seen.insert(key).second) break;

      std::mt19937 gen((uint32_t)(base_seed ^ (0x8F3A2B7C1D9E6A5BULL + i * 1315423911u + tries)));
      Progeny.col(i) = sampleInt_std(nVar, potCross, gen);
      if (++tries > 64) break;
    }
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
  arma::vec x = calcContr(Crosses_mat.rows(idx), nInd);
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

  std::vector<char> keep(N, 1);
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


static inline void cap_archive_by_crowding(
    std::vector<arma::uvec>& plans,
    std::vector<double>& U,
    std::vector<double>& S,
    size_t cap)
{
  const size_t N = U.size();
  if (N <= cap) return;

  std::vector<double> cd(N, 0.0);
  std::vector<size_t> idx_u(N), idx_s(N);
  for (size_t i = 0; i < N; ++i) { idx_u[i] = i; idx_s[i] = i; }

  std::sort(idx_u.begin(), idx_u.end(), [&](size_t a, size_t b){ return U[a] > U[b]; });
  std::sort(idx_s.begin(), idx_s.end(), [&](size_t a, size_t b){ return S[a] < S[b]; });

  double umin = U[idx_u.back()], umax = U[idx_u.front()];
  double smin = S[idx_s.front()], smax = S[idx_s.back()];
  double ur = std::max(umax - umin, 1e-12);
  double sr = std::max(smax - smin, 1e-12);

  cd[idx_u.front()] = cd[idx_u.back()] = std::numeric_limits<double>::infinity();
  cd[idx_s.front()] = cd[idx_s.back()] = std::numeric_limits<double>::infinity();

  for (size_t k = 1; k + 1 < N; ++k) {
    size_t i = idx_u[k];
    if (!std::isinf(cd[i])) cd[i] += (U[idx_u[k-1]] - U[idx_u[k+1]]) / ur;
  }
  for (size_t k = 1; k + 1 < N; ++k) {
    size_t i = idx_s[k];
    if (!std::isinf(cd[i])) cd[i] += (S[idx_s[k+1]] - S[idx_s[k-1]]) / sr;
  }

  std::vector<size_t> order(N);
  for (size_t i = 0; i < N; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b){ return cd[a] > cd[b]; });

  std::vector<char> keep(N, 0);
  for (size_t k = 0; k < cap; ++k) keep[order[k]] = 1;

  std::vector<arma::uvec> plans2;
  std::vector<double> U2, S2;
  plans2.reserve(cap); U2.reserve(cap); S2.reserve(cap);

  for (size_t i = 0; i < N; ++i) if (keep[i]) {
    plans2.push_back(plans[i]);
    U2.push_back(U[i]);
    S2.push_back(S[i]);
  }

  plans.swap(plans2);
  U.swap(U2);
  S.swap(S2);
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
  const arma::uword nVar      = nCross - fixedRows;

  arma::umat Crosses_mat      = as<arma::umat>(Crosses);
  arma::umat fixedCrosses_mat = as<arma::umat>(fixedCrosses);
  Crosses_mat      -= 1;
  fixedCrosses_mat -= 1;
  arma::uword nInd = G.n_cols;


  // contribution of fixed plan
  arma::vec xfixed = calcContr(fixedCrosses_mat, nInd);
  double uMax=0.0, uMin=0.0, simMax=0.0, simMin=0.0;



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

  // --- baseline uMax/simMax using top-u crosses ---
  // Max gain (no GA)
  arma::uvec uMax_plan;
  arma::uvec uBestIndex = arma::sort_index(u, "descend");
  uBestIndex.resize(nVar);
  {
    uMax_plan = uBestIndex; // keep for epsilon init anchor
    arma::vec x = calcContr(Crosses_mat.rows(uBestIndex), nInd);
    uMax  = (arma::accu(u.elem(uBestIndex)) + ufixedSum) / static_cast<double>(nCross);
    simMax = calc_sim(x + xfixed, G);
  }


  // =========================================================
  // Phase 1: GA minimizing sim only
  // =========================================================
  Rcpp::Rcout << "\nPhase 1: Optimize for Similarity (minimize sim)\n";

  ct_parallel_for(0, (int)nPop, [&](int i) {
    std::mt19937 gen((uint32_t)(base_seed ^ (0x9E3779B97F4A7C15ULL + (uint64_t)i)));
    Progeny.col((arma::uword)i) = sampleInt_std(nVar, potCross, gen);

    arma::uvec idx = arma::conv_to<arma::uvec>::from(Progeny.col((arma::uword)i));
    double uu, ss;
    eval_plan(idx, u, ufixedSum, nCross, G, Crosses_mat, nInd, xfixed, uu, ss);

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

    ct_parallel_for(0, (int)nPop, [&](int i) {
      std::mt19937 gen_i((uint32_t)(base_seed ^ (0xD2B74407B1CE6E93ULL + (uint64_t)gen * 1315423911ULL + (uint64_t)i)));

      arma::uvec child = mate_uniform(arma::conv_to<arma::uvec>::from(Parents.col(crossPlan(0, (arma::uword)i))),
                                      arma::conv_to<arma::uvec>::from(Parents.col(crossPlan(1, (arma::uword)i))),
                                      gen_i, potCross);

      std::uniform_real_distribution<> U01(0.0, 1.0);
      if (U01(gen_i) < probMut) child = mutate(child, nMutate, potCross, gen_i);

      Progeny.col((arma::uword)i) = child;

      arma::uvec idx = child;
      double uu, ss;
      eval_plan(child, u, ufixedSum, nCross, G, Crosses_mat, nInd, xfixed, uu, ss);


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


    enforce_unique_population(Progeny, nVar, potCross, base_seed);

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

                                                       ct_parallel_for(0, (int)nPop, [&](int i) {
                                                         std::mt19937 gen_i((uint32_t)(base_seed ^ (0x9E3779B97F4A7C15ULL
                                                                                                       + (uint64_t)gen * 1315423911ULL
                                                                                                       + (uint64_t)i
                                                                                                       + (uint64_t)wi * 0x85EBCA6BULL)));

                                                                                                       arma::uvec child = mate_uniform(
                                                                                                         arma::conv_to<arma::uvec>::from(Parents.col(crossPlan(0, (arma::uword)i))),
                                                                                                         arma::conv_to<arma::uvec>::from(Parents.col(crossPlan(1, (arma::uword)i))),
                                                                                                         gen_i, potCross
                                                                                                       );

                                                                                                       std::uniform_real_distribution<> U01(0.0, 1.0);
                                                                                                       if (U01(gen_i) < probMut) child = mutate(child, nMutate, potCross, gen_i);

                                                                                                       Progeny.col((arma::uword)i) = child;

                                                                                                       double uu, ss;
                                                                                                       eval_plan(child, u, ufixedSum, nCross, G, Crosses_mat, nInd, xfixed, uu, ss);


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

  enforce_unique_population(ProgenyBig, nVar, potCross, base_seed);

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
  enforce_unique_population(Progeny, nVar, potCross, base_seed);

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
  if (ARCH_CAP > 0) cap_archive_by_crowding(arch_plans, arch_u, arch_sim, ARCH_CAP);

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

    ct_parallel_for(0, (int)nPop, [&](int i) {
      std::mt19937 gen_i((uint32_t)(base_seed ^ (0xC949D7C7509E6557ULL + (uint64_t)gen * 40503ULL + (uint64_t)i)));

      arma::uvec child = mate_uniform(arma::conv_to<arma::uvec>::from(Parents.col(crossPlan(0, (arma::uword)i))),
                                      arma::conv_to<arma::uvec>::from(Parents.col(crossPlan(1, (arma::uword)i))),
                                      gen_i, potCross);

      std::uniform_real_distribution<> U01(0.0, 1.0);
      if (U01(gen_i) < probMut) child = mutate(child, nMutate, potCross, gen_i);

      Progeny.col((arma::uword)i) = child;

      double uu, ss;
      eval_plan(child, u, ufixedSum, nCross, G, Crosses_mat, nInd, xfixed, uu, ss);


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
    if (ARCH_CAP > 0) cap_archive_by_crowding(arch_plans, arch_u, arch_sim, ARCH_CAP);

    // parent selection
    std::vector<int> cand2((size_t)nPop);
    std::iota(cand2.begin(), cand2.end(), 0);
    std::vector<int> sel2 = nsga2_select_maxu_minsim(uProgeny, simProgeny, cand2, (int)nSel);

    for (arma::uword i = 0; i < nSel; ++i) {
      int idx = sel2[(size_t)i];
      Parents.col(i) = Progeny.col((arma::uword)idx);
      uParents(i)    = uProgeny((arma::uword)idx);
      simParents(i)  = simProgeny((arma::uword)idx);
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
  arma::umat outCrossPlan = arma::join_cols(Crosses_mat.rows(Best), fixedCrosses_mat);

  // Final Pareto outputs from archive
  Rcpp::DataFrame pareto_df;
  Rcpp::List pareto_crossplans;
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

    for (size_t k = 0; k < nPF; ++k) {
      size_t a = ord[k];
      pu[(int)k] = arch_u[a];
      ps[(int)k] = arch_sim[a];

      arma::uvec vv = arma::sort(arch_plans[a]);
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
    Rcpp::Named("paretoPlans")  = pareto_crossplans
  );
}
