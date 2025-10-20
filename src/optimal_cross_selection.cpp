// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <vector>
#include <cmath>
#include "parallel_backend.h"  // OpenMP when available, else RcppParallel

using namespace Rcpp;

// ------------------------------
// Helpers
// ------------------------------

inline double calc_sim(const arma::vec& x, const arma::mat& G) {
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

// shuffle arma::uvec with given RNG
arma::uvec shuffle_with_rng(const arma::uvec& vec, std::mt19937& gen) {
  std::vector<arma::uword> temp(vec.begin(), vec.end());
  std::shuffle(temp.begin(), temp.end(), gen);
  return arma::uvec(temp);
}

// Map linear index -> row/col in upper triangle of n x n (row <= col)
arma::uword mapRow(const arma::uword& k, const arma::uword& n) {
  return n - 2 - static_cast<arma::uword>(
      std::floor(std::sqrt(-8.0 * static_cast<double>(k)
                             + 4.0 * static_cast<double>(n) * (static_cast<double>(n) - 1.0) - 7.0) / 2.0 - 0.5));
}
arma::uword mapCol(const arma::uword& row, const arma::uword& k, const arma::uword& n) {
  return k + row + 1 - n * (n - 1) / 2 + (n - row) * ((n - row) - 1) / 2;
}

// Thread-safe sampling without replacement: returns sorted unique indices in [0, N-1]
arma::uvec sampleInt_std(arma::uword n, arma::uword N, std::mt19937& gen) {
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
arma::umat sampHalfDialComb_std(arma::uword nLevel, arma::uword n, std::mt19937& gen) {
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
arma::vec calcContr(const arma::umat& crosses, const arma::uword& nInd) {
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

// Normalize to unit square, compute angle and length toward (gain high, similarity low)
void calcVec(double& angle, double& length, double u, double sim,
             const double& uMax, const double& simMax,
             const double& uMin, const double& simMin) {
  u   = (u   - uMin) / (uMax - uMin);
  sim = (simMax - sim) / (simMax - simMin);
  length = std::sqrt(u * u + sim * sim);
  angle  = std::acos(u / std::max(length, 1e-16));
  if (u < 0) length = -length;
}

// Uniform set-based crossover
arma::uvec mate_uniform(const arma::uvec& parent1, const arma::uvec& parent2,
                        std::mt19937& gen, arma::uword potCross) {
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
      if (!present[gene]) { offspring(i) = gene; present[gene] = 1; }
      else offspring(i) = potCross; // mark invalid for now
    }
  }
  for (arma::uword i = 0; i < n; ++i) if (offspring(i) == potCross) {
    arma::uword gene = parent1(i);
    if (!present[gene]) { offspring(i) = gene; present[gene] = 1; }
    else {
      gene = parent2(i);
      if (!present[gene]) { offspring(i) = gene; present[gene] = 1; }
      else {
        for (arma::uword g = 0; g < potCross; ++g) if (!present[g]) { offspring(i) = g; present[g] = 1; break; }
      }
    }
  }
  for (arma::uword i = 0; i < n; ++i) if (offspring(i) >= potCross) {
    for (arma::uword g = 0; g < potCross; ++g) if (!present[g]) { offspring(i) = g; present[g] = 1; break; }
  }
  return offspring;
}

// Mutation via mating with a random plan
arma::uvec mutate(const arma::uvec& crosses, const arma::uword& nMutate,
                  const arma::uword& potCross, std::mt19937& gen) {
  (void)nMutate; // not used in this operator; kept for API parity
  arma::uvec randomPlan = sampleInt_std(crosses.n_elem, potCross, gen);
  return mate_uniform(crosses, randomPlan, gen, potCross);
}

// ------------------------------
// Genetic Algorithm main
// ------------------------------

// [[Rcpp::export]]
Rcpp::List cpp_optimal_cross_selection(const NumericMatrix& Crosses,
                                       const NumericMatrix& fixedCrosses,
                                       arma::uword nCross,
                                       double targetAngle,
                                       arma::vec& u,
                                       arma::vec& ufixed,
                                       arma::mat& G,
                                       double probMut=0.01,
                                       arma::uword nMutate=2,
                                       arma::uword nSel=500,
                                       arma::uword nPop=10000,
                                       arma::uword maxGen=1000,
                                       arma::uword maxRun=100,
                                       double anglePenalty=0.5,
                                       int nThreads=4) {
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

  arma::umat Progeny(nVar, nPop, arma::fill::zeros);
  arma::umat Parents(nVar, nSel, arma::fill::zeros);
  arma::uvec Best(nVar, arma::fill::zeros);

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
  arma::vec xfixed = calcContr(fixedCrosses_mat, nInd);

  // Max gain (no GA)
  arma::uvec uBestIndex = arma::sort_index(u, "descend");
  uBestIndex.resize(nVar);
  {
    arma::vec x = calcContr(Crosses_mat.rows(uBestIndex), nInd);
    uMax  = (arma::accu(u.elem(uBestIndex)) + ufixedSum) / static_cast<double>(nCross);
    simMax = calc_sim(x + xfixed, G);
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

  ct_parallel_for(0, static_cast<int>(nPop), [&](int i) {
    std::mt19937 gen(static_cast<uint32_t>(base_seed ^ (0x9e3779b97f4a7c15ULL + i)));
    Progeny.col(i) = sampleInt_std(nVar, potCross, gen);
    arma::vec x = calcContr(Crosses_mat.rows(Progeny.col(i)), nInd);
    simProgeny(i) = calc_sim(x + xfixed, G);
  });

  rankProgeny = arma::sort_index(simProgeny, "ascend");
  for (arma::uword i = 0; i < nSel; ++i) {
    Parents.col(i) = Progeny.col(rankProgeny(i));
    simParents(i)  = simProgeny(rankProgeny(i));
  }
  simBest = simParents(0);
  Best    = Parents.col(0);

  Rcpp::Rcout << "Gen  Similarity" << std::endl;

  for (arma::uword gen = 0; gen < maxGen; ++gen) {
    arma::umat crossPlan = sampHalfDialComb_std(nSel, nPop, *new std::mt19937(static_cast<uint32_t>(base_seed ^ (0xC3A5C85C97CB3127ULL + gen)))); // local RNG

    ct_parallel_for(0, static_cast<int>(nPop), [&](int i) {
      std::mt19937 gen_i(static_cast<uint32_t>(base_seed ^ (0xD2B74407B1CE6E93ULL + gen * 1315423911u + i)));
      Progeny.col(i) = mate_uniform(Parents.col(crossPlan(0, i)),
                  Parents.col(crossPlan(1, i)), gen_i, potCross);

      std::uniform_real_distribution<> U01(0.0, 1.0);
      if (U01(gen_i) < probMut) {
        Progeny.col(i) = mutate(Progeny.col(i), nMutate, potCross, gen_i);
      }

      arma::vec x = calcContr(Crosses_mat.rows(Progeny.col(i)), nInd);
      simProgeny(i) = calc_sim(x + xfixed, G);
    });

    rankProgeny = arma::sort_index(simProgeny, "ascend");
    for (arma::uword i = 0; i < nSel; ++i) {
      Parents.col(i) = Progeny.col(rankProgeny(i));
      simParents(i)  = simProgeny(rankProgeny(i));
    }

    if (simParents(0) < simBest) {
      simBest = simParents(0);
      Best    = Parents.col(0);
      currentRun = 0;
    } else {
      ++currentRun;
    }

    if (gen % 10 == 0) Rcpp::Rcout << gen << "  " << simBest << std::endl;
    if (currentRun >= maxRun) break;
  }

  simMin  = simBest;
  uMin    = (arma::accu(u.elem(Best)) + ufixedSum) / static_cast<double>(nCross);

  // -------------------------
  // Phase 2: Optimize crossing plan (angle/length w.r.t. targetAngle)
  // -------------------------
  Rcpp::Rcout << std::endl << std::endl << "Optimize for Crossing Plan" << std::endl << std::endl;

  ct_parallel_for(0, static_cast<int>(nPop), [&](int i) {
    std::mt19937 gen_i(static_cast<uint32_t>(base_seed ^ (0xA24BAED4963EE407ULL + i)));
    Progeny.col(i) = sampleInt_std(nVar, potCross, gen_i);

    arma::vec x = calcContr(Crosses_mat.rows(Progeny.col(i)), nInd);
    simProgeny(i) = calc_sim(x + xfixed, G);

    uProgeny(i) = (arma::accu(u.elem(Progeny.col(i))) + ufixedSum) / static_cast<double>(nCross);
    calcVec(angleProgeny(i), lenProgeny(i), uProgeny(i), simProgeny(i),
            uMax, simMax, uMin, simMin);
    valProgeny(i) = lenProgeny(i) - anglePenalty * std::fabs(angleProgeny(i) - targetAngle);
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
    arma::umat crossPlan = sampHalfDialComb_std(nSel, nPop, *new std::mt19937(static_cast<uint32_t>(base_seed ^ (0x9E3779B185EBCA87ULL + gen))));

    ct_parallel_for(0, static_cast<int>(nPop), [&](int i) {
      std::mt19937 gen_i(static_cast<uint32_t>(base_seed ^ (0xC949D7C7509E6557ULL + gen * 40503u + i)));
      Progeny.col(i) = mate_uniform(Parents.col(crossPlan(0, i)),
                  Parents.col(crossPlan(1, i)), gen_i, potCross);

      std::uniform_real_distribution<> U01(0.0, 1.0);
      if (U01(gen_i) < probMut) {
        Progeny.col(i) = mutate(Progeny.col(i), nMutate, potCross, gen_i);
      }

      arma::vec x = calcContr(Crosses_mat.rows(Progeny.col(i)), nInd);
      simProgeny(i) = calc_sim(x + xfixed, G);

      uProgeny(i) = (arma::accu(u.elem(Progeny.col(i))) + ufixedSum) / static_cast<double>(nCross);
      calcVec(angleProgeny(i), lenProgeny(i), uProgeny(i), simProgeny(i),
              uMax, simMax, uMin, simMin);
      valProgeny(i) = lenProgeny(i) - anglePenalty * std::fabs(angleProgeny(i) - targetAngle);
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

  Best = arma::sort(Best);
  arma::umat outCrossPlan = arma::join_cols(Crosses_mat.rows(Best), fixedCrosses_mat);

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
