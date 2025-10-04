// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <omp.h>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <vector>
// #include <chrono>  // Removed as timing is no longer needed

using namespace Rcpp;

//function to shuffle a arma::uvec
arma::uvec shuffle_with_rng(const arma::uvec& vec, std::mt19937& gen) {
    // Convert arma::uvec to std::vector for compatibility with std::shuffle
    std::vector<arma::uword> temp(vec.begin(), vec.end());
    
    // Shuffle using std::shuffle and the provided RNG
    std::shuffle(temp.begin(), temp.end(), gen);
    
    // Convert back to arma::uvec
    return arma::uvec(temp);
}

//Maps a linear index to a row in a symmetric matrix representation.
arma::uword mapRow(const arma::uword& k, const arma::uword& n){
    return n - 2 - static_cast<arma::uword>(std::floor(std::sqrt(-8.0 * static_cast<double>(k) + 4.0 * static_cast<double>(n) * (static_cast<double>(n)-1.0) - 7.0) / 2.0 - 0.5));
}

//Maps a linear index and row to a column in a symmetric matrix representation.
arma::uword mapCol(const arma::uword& row, const arma::uword& k, const arma::uword& n){
    return k + row + 1 - n * (n - 1) / 2 + (n - row) * ((n - row) - 1) / 2;
}

//Randomly samples integers without replacement 
arma::uvec sampleInt(arma::uword n, arma::uword N){
     arma::uvec output;
    output.set_size(n);
    if(n == 0){
        return output;
    }
    double q, v, x, y1, y2;
    arma::uword S, limit, top, bottom;
    arma::vec u(1, arma::fill::randu);
    v = exp(log(u(0)) / static_cast<double>(n));
    q = static_cast<double>(N - n + 1);
    arma::uword threshold = 13 * n;
    while((n > 1) & (threshold < N)){
        while(true){
            while(true){
                x = static_cast<double>(N) * (1.0 - v);
                S = floor(x);
                if(static_cast<double>(S) < q){
                    break;
                }
                u.randu();
                v = exp(log(u(0)) / static_cast<double>(n));
            }
            u.randu();
            y1 = exp(log(u(0) * static_cast<double>(N) / q) / static_cast<double>(n - 1));
            v = y1 * (1.0 - x / static_cast<double>(N)) * (q / (q - static_cast<double>(S)));
            if(v <= 1.0){
                break;
            }
            y2 = 1.0;
            top = N - 1;
            if((n - 1) > S){
                bottom = N - n;
                limit = N - S;
            }
            else{
                bottom = N - S - 1;
                limit = N - n + 1;
            }
            for(arma::uword i = N - 1; i >= limit; --i){
                y2 *= static_cast<double>(top) / static_cast<double>(bottom);
            }
            u.randu();
            if((static_cast<double>(N) / (static_cast<double>(N) - x)) >= (y1 * exp(log(y2) / static_cast<double>(n - 1)))){
                v = exp(log(u(0)) / static_cast<double>(n - 1));
                break;
            }
            v = exp(log(u(0)) / static_cast<double>(n));
        }
        output(n - 1) = S + 1;
        N = N - S - 1;
        --n;
        q = static_cast<double>(N - n + 1);
        threshold -= 13;
    }
    if(n > 1){
        top = N - n;
        while(n >= 2){
            u.randu();
            S = 0;
            q = static_cast<double>(top) / static_cast<double>(N);
            while(q > u(0)){
                ++S;
                --top;
                --N;
                q = (q * static_cast<double>(top)) / static_cast<double>(N);
            }
            output(n - 1) = S + 1;
            --N;
            --n;
        }
        u.randu();
        output(0) = floor(u(0) * static_cast<double>(N));
    }
    else{
        output(0) = floor(v * static_cast<double>(N));
    }
    return arma::cumsum(output);
}

// Samples half diallele combinations without replacement.

arma::umat sampHalfDialComb(arma::uword nLevel, arma::uword n){
    arma::uword N = nLevel * (nLevel - 1) / 2;
    arma::uword fullComb = 0;
    // Determine number of complete combinations
    while(n > N){
        n -= N;
        ++fullComb;
    }
    arma::uvec samples = sampleInt(n, N);
    // Calculate selected combinations
    arma::umat output(2, n);
    for(arma::uword i = 0; i < n; ++i){
        output(0, i) = mapRow(samples(i), nLevel);
        output(1, i) = mapCol(output(0, i), samples(i), nLevel);
    }
    // Add full combinations if necessary
    if(fullComb > 0){
        arma::umat tmp(2, N * fullComb, arma::fill::zeros);
        arma::uword i;
        for(arma::uword j = 0; j < (N * fullComb); ++j){
            i = j % N;
            tmp(0, j) = mapRow(i, nLevel);
            tmp(1, j) = mapCol(tmp(0, j), i, nLevel);
        }
        output = arma::join_cols(output, tmp);
    }
    return output;
}

// Finds parental contributions for a crossing plan.
arma::vec calcContr(const arma::uvec& crosses, const arma::uword& nInd, const arma::uvec& rowPos, const arma::uvec& colPos){
    double val = 1.0 / (2.0 * static_cast<double>(crosses.n_elem));
    arma::vec x(nInd, arma::fill::zeros);
    for(auto it = crosses.begin(); it != crosses.end(); ++it){
        x(rowPos(*it)) += val;
        x(colPos(*it)) += val;
    }
    return x;
}

//Calculates the angle and length of a solution vector, i.e. weights between similarity and gain

void calcVec(double& angle, double& length, double u, double sim,
             const double& uMax, const double& simMax,
             const double& uMin, const double& simMin){
    u = (u - uMin) / (uMax - uMin);
    sim = (simMax - sim) / (simMax - simMin);
    length = std::sqrt(u * u + sim * sim);
    angle = std::acos(u / length);
    if(u < 0){
        // Solution in wrong quadrant
        length = -length;
    }
}

// Performs Uniform Set-Based Crossover between two parents to produce an offspring.
arma::uvec mate_uniform(const arma::uvec& parent1, const arma::uvec& parent2, std::mt19937& gen, arma::uword potCross){
    arma::uword n = parent1.n_elem;
    arma::uvec offspring(n, arma::fill::zeros);
    std::vector<char> present(potCross, 0); // Tracks genes already in offspring
    
    // Initialize offspring with invalid gene markers
    offspring.fill(potCross);
    
    // Uniform crossover: for each gene position, randomly choose from parent1 or parent2
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for(arma::uword i = 0; i < n; ++i){
        bool take_from_p1 = dis(gen) < 0.5;
        arma::uword gene = take_from_p1 ? parent1(i) : parent2(i);
        if(!present[gene]){
            offspring(i) = gene;
            present[gene] = 1;
        }
        else{
            // Try to take the gene from the other parent
            gene = take_from_p1 ? parent2(i) : parent1(i);
            if(!present[gene]){
                offspring(i) = gene;
                present[gene] = 1;
            }
            else{
                // Gene already present; mark as invalid for now
                offspring(i) = potCross;
            }
        }
    }
    
    // Fill any remaining positions with genes from parent1, then parent2
    for(arma::uword i = 0; i < n; ++i){
        if(offspring(i) == potCross){
            // Try to fill from parent1
            arma::uword gene = parent1(i);
            if(!present[gene]){
                offspring(i) = gene;
                present[gene] = 1;
            }
            else{
                // Try to fill from parent2
                gene = parent2(i);
                if(!present[gene]){
                    offspring(i) = gene;
                    present[gene] = 1;
                }
                else{
                    // Find a gene not yet present
                    for(arma::uword g = 0; g < potCross; ++g){
                        if(!present[g]){
                            offspring(i) = g;
                            present[g] = 1;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    // Final check to ensure all genes are unique and within range
    for(arma::uword i = 0; i < n; ++i){
        if(offspring(i) >= potCross){
            // Assign a gene not already present
            for(arma::uword g = 0; g < potCross; ++g){
                if(!present[g]){
                    offspring(i) = g;
                    present[g] = 1;
                    break;
                }
            }
        }
    }
    
    return offspring;
}

// Performs mutation of a plan by mating it with a random crossing plan.
 
 arma::uvec mutate(const arma::uvec& crosses, const arma::uword& nMutate, const arma::uword& potCross, std::mt19937& gen){
    arma::uvec mutations = sampleInt(nMutate, potCross);
    // Create a random crossing plan for mating
    arma::uvec randomPlan = sampleInt(crosses.n_elem, potCross);
    arma::uvec mutated = mate_uniform(crosses, randomPlan, gen, potCross);
    
    return mutated;
}

//Genetic Algorithm to select optimal crossing plans.
// [[Rcpp::export]]
Rcpp::List cpp_optimal_cross_selection(arma::uword nCross, // Number of crosses to make
                         double targetAngle, // Target angle betwen maximum gain and diversity (in radians)
                         arma::vec& u, // Vector of criteria for crosses
                         arma::mat& G, // Relationship matrix among individuals 
                         double probMut=0.01, // Mutation probability for progeny in GA
                         arma::uword nMutate=2, // Number of potential mutations in mutated progeny
                         arma::uword nSel=500, // Number of parents in GA
                         arma::uword nPop=10000, // Number of progeny in GA
                         arma::uword maxGen=1000, // Maximum number of generations
                         arma::uword maxRun=100, // Stopping criteria for maximum number of runs without change
                         double anglePenalty=0.5, // Penalty to vector length for off angle, higher value emphasizes angle more
                         int nThreads=4){ // Number of threads for OpenMP
    omp_set_num_threads(nThreads); // Sets number of threads for OpenMP
    
    // Initialize rowPos and colPos mappings
    arma::uword nInd = G.n_cols; // Number of individuals
    arma::uword potCross = nInd * (nInd - 1) / 2; // Potential crosses
    arma::uvec rowPos(potCross, arma::fill::zeros), colPos(potCross, arma::fill::zeros);
    arma::uword tmpI = 0;
    for(arma::uword i = 0; i < (nInd - 1); i++){
        for(arma::uword j = (i + 1); j < nInd; j++){
            rowPos(tmpI) = i;
            colPos(tmpI) = j;
            ++tmpI;
        }
    }
    
    // Initialize progeny and parents matrices
    arma::umat Progeny(nCross, nPop, arma::fill::zeros), Parents(nCross, nSel, arma::fill::zeros);
    arma::uvec Best(nCross, arma::fill::zeros);
    arma::vec uProgeny(nPop, arma::fill::zeros), uParents(nSel, arma::fill::zeros);
    arma::vec simProgeny(nPop, arma::fill::zeros), simParents(nSel, arma::fill::zeros);
    arma::vec angleProgeny(nPop, arma::fill::zeros), angleParents(nSel, arma::fill::zeros);
    arma::vec lenProgeny(nPop, arma::fill::zeros), lenParents(nSel, arma::fill::zeros);
    arma::vec valProgeny(nPop, arma::fill::zeros), valParents(nSel, arma::fill::zeros);
    arma::uvec rankProgeny(nPop, arma::fill::zeros);
    double uBest_val, simBest, valBest, angleBest, lenBest;
    double uMax, uMin, simMax, simMin;
    arma::uword currentRun = 0;
    arma::umat crossPlan(2, nPop, arma::fill::zeros); // Declare crossPlan here
    
    // Initialize random number generator for initial population
    std::random_device rd_initial;
    std::mt19937 gen_initial(rd_initial());
    
    // Calculate maximum gain (doesn't require GA)
    arma::uvec uBestIndex = arma::sort_index(u, "descend");
    uBestIndex.resize(nCross);
    arma::vec x = calcContr(uBestIndex, nInd, rowPos, colPos);
    uMax = arma::mean(u.elem(uBestIndex));
    simMax = arma::as_scalar(x.t() * G * x);
    
    if(targetAngle < 1e-6){
        // No need to optimize crossing plan, just return the best
        arma::umat outCrossPlan(nCross, 2, arma::fill::zeros);
        for(arma::uword i = 0; i < nCross; ++i){
            outCrossPlan(i, 0) = rowPos(uBestIndex(i));
            outCrossPlan(i, 1) = colPos(uBestIndex(i));
        }
        return Rcpp::List::create(
            Rcpp::Named("crossPlan") = outCrossPlan + 1, // Adjust for 1-based indexing in R
            Rcpp::Named("uMax") = uMax,
            Rcpp::Named("simMax") = simMax
            // Timing information removed
        );
    }
    
    // Optimize for Minimum Similarity
    Rcpp::Rcout << "Optimize for Similarity" << std::endl << std::endl;
    
    // Initialize progeny in parallel with thread-local RNGs
    #pragma omp parallel
    {
        // Each thread gets its own RNG instance
        std::mt19937 gen_thread(rd_initial() + omp_get_thread_num());
        
        #pragma omp for schedule(static)
        for(arma::uword i = 0; i < nPop; i++){
            Progeny.col(i) = sampleInt(nCross, potCross);
            arma::vec x = calcContr(Progeny.col(i), nInd, rowPos, colPos);
            simProgeny(i) = arma::as_scalar(x.t() * G * x);
        }
    }
    
    // Select parents and best solution based on similarity
    rankProgeny = arma::sort_index(simProgeny, "ascend");
    for(arma::uword i = 0; i < nSel; i++){
        Parents.col(i) = Progeny.col(rankProgeny(i));
        simParents(i) = simProgeny(rankProgeny(i));
    }
    simBest = simParents(0);
    Best = Parents.col(0);
    
    // Run GA for similarity optimization
    Rcpp::Rcout << "Gen  Similarity" << std::endl;
    
    for(arma::uword gen = 0; gen < maxGen; gen++){
        // Generate mating pairs
        crossPlan = sampHalfDialComb(nSel, nPop);
        
        // Parallel region with thread-local RNGs
        #pragma omp parallel
        {
            // Each thread gets its own RNG instance
            std::mt19937 gen_thread(rd_initial() + omp_get_thread_num());
                
            #pragma omp for schedule(static)
            for(arma::uword i = 0; i < nPop; i++){
                Progeny.col(i) = mate_uniform(Parents.col(crossPlan(0, i)),
                                             Parents.col(crossPlan(1, i)), gen_thread, potCross);
                arma::vec p(1, arma::fill::randu);
                if(as_scalar(p) < probMut){
                    Progeny.col(i) = mutate(Progeny.col(i), nMutate, potCross, gen_thread);
                }
                arma::vec x = calcContr(Progeny.col(i), nInd, rowPos, colPos);
                simProgeny(i) = arma::as_scalar(x.t() * G * x);
            }
        }
        
        // Select the best parents based on updated similarity
        rankProgeny = arma::sort_index(simProgeny, "ascend");
        for(arma::uword i = 0; i < nSel; i++){
            Parents.col(i) = Progeny.col(rankProgeny(i));
            simParents(i) = simProgeny(rankProgeny(i));
        }
        
        // Update the best similarity if improved
        if(simParents(0) < simBest){
            simBest = simParents(0);
            Best = Parents.col(0);
            currentRun = 0;
        }
        else{
            ++currentRun;
        }
        
        // Report status every 10 generations
        if(gen % 10 == 0){
            Rcpp::Rcout << gen << "  " << simBest << std::endl;
        }
        
        // Terminate if no improvement for maxRun generations
        if(currentRun >= maxRun){
            break;
        }
    }
    
    // Record minimum similarity and associated gain
    simMin = simBest;
    uMin = arma::mean(u.elem(Best));
    
    // Optimize for Crossing Plan
    Rcpp::Rcout << std::endl << std::endl << "Optimize for Crossing Plan" << std::endl << std::endl;
    
    // Initialize progeny for crossing plan optimization in parallel with thread-local RNGs
    #pragma omp parallel
    {
        // Each thread gets its own RNG instance
        std::mt19937 gen_thread(rd_initial() + omp_get_thread_num());
            
        #pragma omp for schedule(static)
        for(arma::uword i = 0; i < nPop; i++){
            Progeny.col(i) = sampleInt(nCross, potCross);
            arma::vec x = calcContr(Progeny.col(i), nInd, rowPos, colPos);
            simProgeny(i) = arma::as_scalar(x.t() * G * x);
            uProgeny(i) = arma::mean(u.elem(Progeny.col(i)));
            calcVec(angleProgeny(i), lenProgeny(i), uProgeny(i),
                    simProgeny(i), uMax, simMax, uMin, simMin);
            valProgeny(i) = lenProgeny(i) - anglePenalty * std::fabs(angleProgeny(i) - targetAngle);
        }
    }
    
    // Select parents based on crossing plan value
    rankProgeny = arma::sort_index(valProgeny, "descend");
    for(arma::uword i = 0; i < nSel; i++){
        Parents.col(i) = Progeny.col(rankProgeny(i));
        uParents(i) = uProgeny(rankProgeny(i));
        simParents(i) = simProgeny(rankProgeny(i));
        angleParents(i) = angleProgeny(rankProgeny(i));
        lenParents(i) = lenProgeny(rankProgeny(i));
        valParents(i) = valProgeny(rankProgeny(i));
    }
    
    valBest = valParents(0);
    angleBest = angleParents(0);
    lenBest = lenParents(0);
    uBest_val = uParents(0);
    simBest = simParents(0);
    Best = Parents.col(0);
    
    // Run GA for crossing plan optimization
    Rcpp::Rcout << "Gen  acc_cross_valie  Similarity  Angle  Length  Value" << std::endl;
    currentRun = 0;
    
    for(arma::uword gen = 0; gen < maxGen; gen++){
        // Generate mating pairs
        crossPlan = sampHalfDialComb(nSel, nPop);
        
        // Parallel region with thread-local RNGs
        #pragma omp parallel
        {
            // Each thread gets its own RNG instance
            std::mt19937 gen_thread(rd_initial() + omp_get_thread_num());
                
            #pragma omp for schedule(static)
            for(arma::uword i = 0; i < nPop; i++){
                Progeny.col(i) = mate_uniform(Parents.col(crossPlan(0, i)),
                                             Parents.col(crossPlan(1, i)), gen_thread, potCross);
                arma::vec p(1, arma::fill::randu);
                if(as_scalar(p) < probMut){
                    Progeny.col(i) = mutate(Progeny.col(i), nMutate, potCross, gen_thread);
                }
                arma::vec x = calcContr(Progeny.col(i), nInd, rowPos, colPos);
                simProgeny(i) = arma::as_scalar(x.t() * G * x);
                uProgeny(i) = arma::mean(u.elem(Progeny.col(i)));
                calcVec(angleProgeny(i), lenProgeny(i), uProgeny(i),
                        simProgeny(i), uMax, simMax, uMin, simMin);
                valProgeny(i) = lenProgeny(i) - anglePenalty * std::fabs(angleProgeny(i) - targetAngle);
            }
        }
        
        // Select the best parents based on updated crossing plan value
        rankProgeny = arma::sort_index(valProgeny, "descend");
        for(arma::uword i = 0; i < nSel; i++){
            Parents.col(i) = Progeny.col(rankProgeny(i));
            uParents(i) = uProgeny(rankProgeny(i));
            simParents(i) = simProgeny(rankProgeny(i));
            angleParents(i) = angleProgeny(rankProgeny(i));
            lenParents(i) = lenProgeny(rankProgeny(i));
            valParents(i) = valProgeny(rankProgeny(i));
        }
        
        // Update the best value if improved
        if(valParents(0) > valBest){
            valBest = valParents(0);
            angleBest = angleParents(0);
            lenBest = lenParents(0);
            uBest_val = uParents(0);
            simBest = simParents(0);
            Best = Parents.col(0);
            currentRun = 0;
        }
        else{
            ++currentRun;
        }
        
        // Report status every 10 generations
        if(gen % 10 == 0){
            Rcpp::Rcout << gen << "  " << uBest_val << "  " << simBest << "  " 
                       << angleBest << "  " << lenBest << "  " << valBest << std::endl;
        }
        
        // Terminate if no improvement for maxRun generations
        if(currentRun >= maxRun){
            break;
        }
    }
    
    // Convert the best solution to an ordered crossing plan
    Best = arma::sort(Best);
    arma::umat outCrossPlan(nCross, 2, arma::fill::zeros);
    for(arma::uword i = 0; i < nCross; ++i){
        outCrossPlan(i, 0) = rowPos(Best(i));
        outCrossPlan(i, 1) = colPos(Best(i));
    }
    
    // Return the results as an R list, excluding timing information
    return Rcpp::List::create(
        Rcpp::Named("crossPlan") = outCrossPlan + 1, // Adjust for 1-based indexing in R
        Rcpp::Named("uMax") = uMax,
        Rcpp::Named("uMin") = uMin,
        Rcpp::Named("simMax") = simMax,
        Rcpp::Named("simMin") = simMin,
        Rcpp::Named("uBest") = uBest_val,
        Rcpp::Named("simBest") = simBest,
        Rcpp::Named("angleBest") = angleBest,
        Rcpp::Named("lenBest") = lenBest
        // Timing information removed
    );
}
