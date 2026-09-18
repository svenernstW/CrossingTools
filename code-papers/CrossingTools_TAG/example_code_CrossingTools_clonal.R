###############################################################################
# Libraries
###############################################################################

library(AlphaSimR)
library(asreml)
library(CrossingTools)

set.seed(123)

###############################################################################
# Create founders
###############################################################################

###############################################################################
# Create founders
###############################################################################

# 10 chromosomes
# 100 Mb per chromosome
# 100 cM = 1 Morgan per chromosome
# 50 diploid, outbred founders
# 1000 segregating sites retained per chromosome:
#   500 will become QTL
#   500 will become SNP markers

founderPop <- runMacs(
  nInd     = 50,
  nChr     = 10,
  segSites = 1000,
  inbred   = FALSE,
  ploidy   = 2
)

###############################################################################
# Initialize simulation parameters
###############################################################################

SP <- SimParam$new(founderPop)

###############################################################################
# SNP chip
###############################################################################

SP$addSnpChip(
  nSnpPerChr = 500
)

###############################################################################
# Trait
###############################################################################

#loosly following Werner et al. (2023) (https://doi.org/10.1007/s00122-023-04300-6)

SP$addTraitAD(
  nQtlPerChr = 500,     
  mean        = 0,
  var         = 1,       
  meanDD      = 0.3,
  varDD       = 0.2,
  useVarA     = TRUE
)


###############################################################################
# Create founder parents and assign phenotypes (EYT evaluation)
###############################################################################

# Create founder parents
Parents = newPop(founderPop)

# Add error to get a phenotype with desired H2
Parents = setPheno(Parents, h2 = 0.5)

rm(founderPop)

###############################################################################
# Create outbred training population
###############################################################################

F1 <- randCross(
  Parents,
  nCrosses = 50,
  nProgeny = 20,
  simParam = SP
)

pop <- setPheno(
  F1,
  H2 = 0.5,
  simParam = SP
)

###############################################################################
# Prepare genomic and phenotypic data for additive-dominance GBLUP
###############################################################################

# Marker genotype matrix for the training population
Mtrain <- pullSnpGeno(pop)

# Haplotype matrix, in reality, Genotype data would need to be phased with tools such as BEAGLE, only needed for SPV and segregation variances
Hap1 <- pullSnpHaplo(pop)[seq(1,nrow(Mtrain)*2,2),]
Hap2 <- pullSnpHaplo(pop)[seq(2,nrow(Mtrain)*2,2),]

row.names(Hap1) <- row.names(Mtrain)
row.names(Hap2) <- row.names(Mtrain)

# Phenotype table for genomic prediction
Pheno <- data.frame(
  id    = as.factor(pop@id),
  value = pop@pheno[,1]
)

###############################################################################
# Additive genomic relationship matrix
###############################################################################

# Allele frequencies
p <- colMeans(Mtrain) / 2
q <- 1 - p

# Additive coding: W = M - 2p
W <- Mtrain - 2 * p

# Additive GRM
GRM.A <- tcrossprod(W) / ncol(Mtrain)


###############################################################################
# Dominance-deviation genomic relationship matrix
###############################################################################

# Heterozygosity indicator: 1 for heterozygotes, 0 otherwise
H <- 1 * (Mtrain == 1)

# Statistical dominance-deviation coding:
# Z = H - (1 - 2p)W - 2pq
Z <- H - (1 - 2 * p) * W - 2 * p * q

# Dominance GRM
GRM.D <- tcrossprod(Z) / ncol(Mtrain)

###############################################################################
# Fit additive-dominance GBLUP model (ASReml)
###############################################################################

# Increase ASReml workspace
asreml.options(pworkspace = "5gb",workspace = "5gb")

GBLUP <- asreml(fixed = value ~ 1,
                random = ~vm(id,GRM.A)+vm(id,GRM.D),
                residual = ~ id(units),
                na.action = na.method(y = "include"),
                data = Pheno)





###############################################################################
# Extract breeding values and dominance deviations
###############################################################################

# Breeding values for Trait1 and Trait2 (stacked random coefficients)
A <- as.data.frame(GBLUP$coefficients$random[1:1000,])
D <- as.data.frame(GBLUP$coefficients$random[1001:2000,])


names(A) <- c("Trait1")
names(D) <- c("Trait1")

###############################################################################
# Back-solve additive and dominance marker effects
###############################################################################

# Convert breeding values to marker effects via backsolve_marker_effects()
alpha <- backsolve_marker_effects(
  marker.mat     = W,
  G.mat          = GRM.A,
  genotype.effects        = A
)

head(alpha)

# Convert dominance deviation to marker effects via backsolve_marker_effects()
delta <- backsolve_marker_effects(
  marker.mat     = Z,
  G.mat          = GRM.D,
  genotype.effects        = D
)

head(delta)
###############################################################################
# Genetic map for markers
###############################################################################

# Get marker genetic map and assign sequential site IDs
map        <- getSnpMap()
map$site   <- seq_len(nrow(map))
#map$site   <- map$id

###############################################################################
# Strategy:  cross expectation (average GEBV and predicted F1 performance), Optimal haploid value and superior progeny index (including additive variance and dominance)
###############################################################################

# Build all possible parent crosses
PotCrosses <- make_cross_plan(parents = 1:nrow(Mtrain))
head(PotCrosses)

# Expected breeding value and total genotypic value of each F1 cross
expectations <- calc_midparent_outcross(
  crosses         = PotCrosses,
  marker.mat      = Mtrain,
  marker.effects.A         = alpha,
  marker.effects.D         = delta,
  nthreads       = 7
)

head(expectations)

# Optimal haploid value for each cross
ohv <- calc_optimal_haploid_value(
  crosses         = PotCrosses,
  marker.mat      = Mtrain,
  marker.effects         = alpha,
  nthreads       = 7
)

head(ohv)

# Additive and dominance segregation variance and superior progeny value
# Selection intensity at cross level
nCrosses <- 50
a     <- nCrosses / nrow(PotCrosses) # this is really just a weighting factor on the segregation variance
intensity <- dnorm(qnorm(1 - a), 0, 1) / a


#with two traits and about two million crosses as in this example might take a while
spv <- calc_spv_outcross(
  crosses         = PotCrosses,
  genetic.map     = map,
  hap.mat1        = Hap1,
  hap.mat2        = Hap2,
  marker.effects.A= alpha,
  marker.effects.D= delta,
  intensity       = intensity,
  covariance      = F,
  nthreads        = 10)

head(spv)


###############################################################################
# Optimal cross selection (trade-off: diversity vs gain; here using SPV of index)
###############################################################################

# Notes:
# - Can add fixed.crosses (always conducted)
# - Can remove potential crosses from optimization by removing them from `crosses`
ocs_pareto <- optimize_cross_plan(
  candidate.crosses       = PotCrosses,
  criterion             = spv$TSPV.Trait1,
  G.mat             = GRM.A,
  method = "pareto", #either pareto to return a pareto with multiple solutions or "angle" to only find a single solution to maximize gain and balance diversity alon a given target angle
  ncrosses      = 50,
  plot=T
)

ocs_pareto$pareto.plans[[1]] # a list of all pareto plans
head(ocs_pareto$pareto.frontier) # coordinates on the pareto

#alternativly, optimize a long a single target angle between gain and genetic similarity

ocs_angle <- optimize_cross_plan(
  candidate.crosses       = PotCrosses,
  criterion             = spv$TSPV.Trait1,
  G.mat             = GRM.A,
  method = "angle", #either pareto to return a pareto with multiple solutions or "angle" to only find a single solution to maximize gain and balance diversity alon a given target angle
  ncrosses      = 50,
  target.angle  = 15
)
head(ocs_angle)





###############################################################################
# Final cross plan and evaluation
###############################################################################

# Final cross plan
crosses <- ocs_pareto$pareto.plans[[3000]] # pick the one that fits your goals best


###############################################################################
# Evaluate and plot a crossing plan
###############################################################################
# Evaluate the plan

cross.df <- cbind(spv,ohv[,3])


summarize_cross_plan(
  cross.plan  = crosses,
  cross.df = cross.df
)

# Plot the plan
plot_cross_plan(
  cross.plan  = crosses,
  cross.df = cross.df
)
