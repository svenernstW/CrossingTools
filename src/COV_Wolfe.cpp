// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <thread> // For sleep functionality
#include <chrono> // For time units

using namespace Rcpp;
using namespace arma;

// Inline function to calculate recombination fraction
inline double qjk(double x, double y) {
  double diff = std::abs(x - y);
  double rcf = 0.5 * (1 - std::exp(-2 * diff));
  return (1 - 2 * rcf);
}

// [[Rcpp::export]]
SEXP cpp_calculate_covariance_wolfe(const NumericMatrix& Crosses,
             const List& genMap,
             const NumericMatrix& Hap1, 
             const NumericMatrix& Hap2,
             const NumericMatrix& U,
             const NumericMatrix& D,
             double intensity,
             bool covariance = false,
             int nThreads = 4) {
#ifdef _OPENMP
  omp_set_num_threads(nThreads);  // set OpenMP threads
#endif
  
  const arma::uword numCrosses   = Crosses.nrow();
  const arma::uword numMarkers   = Hap1.ncol();
  const arma::uword numTrait     = U.ncol();
  const arma::uword numTraitComb = numTrait * (numTrait + 1) / 2;
  
  // Convert inputs to Armadillo
  arma::mat Hap1_mat = as<arma::mat>(Hap1);   // (n_individuals × numMarkers)
  arma::mat Hap2_mat = as<arma::mat>(Hap2);   // (n_individuals × numMarkers)
  arma::mat M_mat = Hap1_mat+Hap2_mat;
  arma::mat U_mat = as<arma::mat>(U);   // (numMarkers × numTrait)
  arma::mat D_mat = as<arma::mat>(D);   // (numMarkers × numTrait)
  
  const arma::uword OFF_EG  = 0;
  const arma::uword OFF_ETG  = numTrait;
  const arma::uword OFF_VAR_A = 2 * numTrait;
  const arma::uword OFF_SPV = 3 * numTrait;
  const arma::uword OFF_VAR_D = 4 * numTrait;
  const arma::uword OFF_TSPV = 5 * numTrait;
  
  
  // Precompute chromosome ranges
  
  std::vector<std::pair<int, int>> chromosomeRanges;
  {
    int startIdx = 0;
    for (int i = 0; i < genMap.size(); ++i) {
      
      NumericMatrix genMapMatrix = as<NumericMatrix>(genMap[i]);
      int n = genMapMatrix.nrow();
      chromosomeRanges.emplace_back(startIdx, startIdx + n - 1);
      startIdx += n;
    }
  }
  
  // Precompute recombination fractions for each generation scenario
  
  mat RC(numMarkers, numMarkers, fill::zeros);
  {
    int startIdx = 0;
    for (int i = 0; i < genMap.size(); ++i) {
      
      NumericMatrix genMapMatrix = as<NumericMatrix>(genMap[i]);
      int n = genMapMatrix.nrow();
      for (int j = 0; j < n; ++j) {
        for (int k = j; k < n; ++k) { 
          double value = qjk(genMapMatrix(j, 0), genMapMatrix(k, 0));
          RC(startIdx + j, startIdx + k) = value;
          RC(startIdx + k, startIdx + j) = value; 
        }
      }
      startIdx += n;
    }
  }
  
  
  
  
  // Copy parent indices out of Crosses (1-based in R → 0-based here)
  arma::Col<int> P1_idx(numCrosses), P2_idx(numCrosses);
  for (arma::uword x = 0; x < numCrosses; ++x) {
    P1_idx[x] = static_cast<int>(Crosses(x, 0)) - 1;
    P2_idx[x] = static_cast<int>(Crosses(x, 1)) - 1;
  }
  
  // Results
  arma::mat results1A(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results1D(numCrosses, numTraitComb, arma::fill::zeros);
  arma::mat results2(numCrosses, numTrait*6, arma::fill::zeros);
  
  // Helper to map (ti, tj) with 0 ≤ ti ≤ tj < numTrait to column index
  auto tri_u_idx_incl = [numTrait](arma::uword ti, arma::uword tj) -> arma::uword {
    return ti * numTrait - (ti * (ti - 1)) / 2 + (tj - ti);
  };
  
#pragma omp parallel for schedule(dynamic)
  for (R_xlen_t x = 0; x < (R_xlen_t)numCrosses; ++x) {
    const int P1 = P1_idx[x];
    const int P2 = P2_idx[x];
    
    arma::uvec differingMarkers = arma::find(
      (M_mat.row(P1) == 1.0) + (M_mat.row(P2) == 1.0)
    );
    
    for (arma::uword ti = 0; ti < numTrait; ++ti) {
      for (arma::uword tj = ti; tj < numTrait; ++tj) {
        if(!covariance && ti!= tj) continue;
        
        double SigmaSqAP1P2 = 0.0;
        double SigmaSqDP1P2 = 0.0;
        
        for (const auto& chrRange : chromosomeRanges) {
          int startC = chrRange.first;
          int endC = chrRange.second;
          
          arma::uvec chrDiff = differingMarkers( arma::find(
            (differingMarkers >= (arma::uword)startC) % (differingMarkers <= (arma::uword)endC)  // element-wise AND
          ));
          
          if (chrDiff.n_elem == 0) {
            continue;
          }
          
          
          for (uword i = 0; i < chrDiff.n_elem; ++i) {
            int marker_i = (int)chrDiff[i];
            for (uword j = i; j < chrDiff.n_elem; ++j) {
              int marker_j = (int)chrDiff[j];
              // Subset recombination fractions
              double rcf_ij = RC(marker_i, marker_j);
              
              arma::uword row_idx_P1  = static_cast<arma::uword>(P1);
              arma::uword row_idx_P2  = static_cast<arma::uword>(P2);
              
              arma::uword col_idx1 = static_cast<arma::uword>(marker_i);
              arma::uword col_idx2 = static_cast<arma::uword>(marker_j);
              
              arma::uvec cols(2);
              cols(0) = col_idx1;
              cols(1) = col_idx2;
              
              // Extract rows for marker_i and marker_j
              arma::rowvec Hap1_P1 = Hap1_mat(arma::uvec{row_idx_P1}, cols);
              arma::rowvec Hap2_P1 = Hap2_mat(arma::uvec{row_idx_P1}, cols);
              
              arma::rowvec Hap1_P2 = Hap1_mat(arma::uvec{row_idx_P2}, cols);
              arma::rowvec Hap2_P2 = Hap2_mat(arma::uvec{row_idx_P2}, cols);
              
              arma::rowvec P1_means = 0.5 * (Hap1_P1 + Hap2_P1);
              arma::rowvec P2_means = 0.5 * (Hap1_P2 + Hap2_P2);
              
              arma::mat Hap_P1 = join_cols(Hap1_P1, Hap2_P1);
              arma::mat Hap_P2 = join_cols(Hap1_P2, Hap2_P2);
              
              mat D1 = 0.5 * Hap_P1.t() * Hap_P1 - P1_means.t() * P1_means;
              mat D2 = 0.5 * Hap_P2.t() * Hap_P2 - P2_means.t() * P2_means;
              
              
              //Compute DGen for the pair (marker_i, marker_j)
              double DGen_ij = (rcf_ij * D1.at(0,1)) + (rcf_ij *D2.at(0,1));
              
              
              
              double contribA = U_mat(marker_j,tj) * DGen_ij * U_mat(marker_i,ti);
              SigmaSqAP1P2 += (i == j) ? contribA : 2 * contribA;
              
              
              double contribD = D_mat(marker_j,tj) * (DGen_ij * DGen_ij) * D_mat(marker_i,ti);
              SigmaSqDP1P2 += (i == j) ? contribD : 2 * contribD;
              
              
              
            }
          }
        }
        
        const arma::uword k = tri_u_idx_incl(ti, tj);
        results1A(x, k) = SigmaSqAP1P2;
        results1D(x, k) = SigmaSqDP1P2;
        
        if (ti == tj) {
          double G1 = arma::dot(M_mat.row(P1), U_mat.col(ti));
          double G2 = arma::dot(M_mat.row(P2), U_mat.col(ti));
          
          double eG = 0.5 * (G1 + G2);
           
          arma::rowvec pk = M_mat.row(P1);      
          arma::rowvec qk = 2.0 - pk;                 
          arma::rowvec ykl = pk - M_mat.row(P2);
          const arma::colvec v1 = (pk - qk - ykl).t();
          const arma::colvec v2 = (2.0 * (pk % qk) + (ykl % (pk - qk))).t();
          
          double eTG = arma::dot(U_mat.col(ti), v1)
            + arma::dot(D_mat.col(ti), v2);
          
          
          results2(x, OFF_EG  + ti) = eG;
          results2(x, OFF_VAR_A + ti) = SigmaSqAP1P2;
          results2(x, OFF_SPV + ti) = eG + intensity * std::sqrt(SigmaSqAP1P2);
          results2(x, OFF_VAR_D + ti) = SigmaSqDP1P2;
          results2(x, OFF_ETG + ti) = eTG;
          results2(x, OFF_TSPV + ti) = eTG + (std::sqrt(SigmaSqAP1P2)+std::sqrt(SigmaSqDP1P2));
          
          
        }
      }
    }
  }
  // If requested, convert each row to an nTrait×nTrait symmetric matrix
  if (covariance) {
    List out_A(numCrosses);
    List out_D(numCrosses);
    for (uword x = 0; x < numCrosses; ++x) {
      NumericMatrix GA(numTrait, numTrait);
      NumericMatrix GD(numTrait, numTrait);
      for (uword ti = 0; ti < numTrait; ++ti) {
        for (uword tj = ti; tj < numTrait; ++tj) {
          const uword k = tri_u_idx_incl(ti, tj);
          double va = results1A(x, k);
          GA(ti, tj) = va;
          GA(tj, ti) = va;  // mirror
          double vd = results1D(x, k);
          GD(ti, tj) = vd;
          GD(tj, ti) = vd;  // mirror
        }
      }
      out_A[x] = GA;
      out_D[x] = GD;
    }
    return List::create(
      Named("cross_values") = results2,
      Named("covariances") = List::create(Named("additive") = out_A, Named("dominance") = out_D)
    );
  }
  
  return Rcpp::wrap(results2);    // returns upper triangle of the covariance matrix in vectorized form
}
