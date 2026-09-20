###############################################################################
# Libraries
###############################################################################
library(AlphaSimR)
library(asreml)
library(CrossingTools)
###############################################################################
# Create founders
###############################################################################

# Generate initial haplotypes
founderPop = runMacs(nInd     = 50,
                     nChr     = 10,
                     segSites = 1000,
                     inbred   = TRUE,
                     species  = "WHEAT")

# Initialize simulation parameters from founder population
SP = SimParam$new(founderPop)

# Define SNP chip density
SP$addSnpChip(500)

###############################################################################
# Add traits
###############################################################################

# Add (additive) multi-trait architecture (e.g., yield-related traits)
SP$addTraitA(nQtlPerChr = 500,
             mean       = c(0,0),
             var        = c(1,2),
             corA       = matrix(c(1,-0.5,
                                   -0.5,1),nrow=2))


###############################################################################
# Create founder parents and assign phenotypes (EYT evaluation)
###############################################################################

# Create founder parents
Parents = newPop(founderPop)

# Add error to get a phenotype with desired H2
Parents = setPheno(Parents, H2 = c(0.5,0.5))

rm(founderPop)

###############################################################################
# Create F1 families and doubled haploid training population
###############################################################################

F1 <- randCross(Parents, 50)
DH <- makeDH(F1, 40)

# Phenotype DH lines
DH <- setPheno(DH, H2 = c(0.5, 0.5))
###############################################################################
# Prepare genomic + phenotypic data for multi-trait GBLUP
###############################################################################

# Marker matrix for training (DH)
Mtrain <- pullSnpGeno(DH)

# Long-format phenotype table for multi-trait model
Pheno <- data.frame(
  id    = as.factor(rep(DH@id,2)),
  trait = as.factor(c(rep("Trait1",nInd(DH)),rep("Trait2",nInd(DH)))),
  value = c(DH@pheno[,1],DH@pheno[,2])
)

###############################################################################
# Additive genomic relationship matrix
###############################################################################

# Allele frequencies used to centre marker genotypes
p <- colMeans(Mtrain) / 2

# Additive genotype coding: W = M - 2p
W <- Mtrain - 2 * p

# Additive genomic relationship matrix
GRM <- tcrossprod(W) / ncol(Mtrain)

###############################################################################
# Fit GBLUP model (ASReml) and predict breeding values
###############################################################################

# Increase ASReml workspace
asreml.options(pworkspace = "4gb",workspace = "4gb")

GBLUP <- asreml(fixed = value ~ 1+trait,
                random = ~vm(id,GRM):us(trait),
                residual = ~ dsum(~ id(units) | trait),
                na.action = na.method(y = "include"),
                data = Pheno)




# EBVs (BLUPs) for id x trait, plus vcov matrix of those predictions

pred <- predict(
  GBLUP,
  classify = "vm(id, GRM):trait",
  only = "vm(id, GRM):trait",
  vcov = TRUE
)


# Prediction-error covariance matrix of the breeding-value BLUPs
PEV <- as.matrix(pred$vcov)

G <- diag(GBLUP$vparameters[c(1,3)],2,2) #this is the genetic covariance matrix
G[2,1] <- G[1,2] <- GBLUP$vparameters[c(2)]


# Prior covariance matrix of the stacked breeding values
GRM.G <- GRM %x% G

# Covariance matrix of the predicted breeding values:
# Var(A_hat) = Var(A) - PEV
V <- GRM.G - PEV

###############################################################################
# Extract breeding values (two traits)
###############################################################################

# Breeding values for Trait1 and Trait2 (stacked random coefficients)
A <- as.data.frame(cbind(
  GBLUP$coefficients$random[grepl("Trait1",row.names(GBLUP$coefficients$random))],
  GBLUP$coefficients$random[grepl("Trait2",row.names(GBLUP$coefficients$random))]
))

names(A) <- c("Trait1","Trait2")

###############################################################################
# Summarize expected index behavior from trait covariance (quick intuition)
###############################################################################
desired_gain <- c(10, 10)  # desired relative gain of 1:1 across traits
smith_hazel  <- c(1, 1)    # equal economic weights for both traits

# Uses var.mat (trait covariance matrix), not per-genotype VCOV
predict_response(var.mat = cov(A), desired.gains = desired_gain, intensity = 1)

predict_response(var.mat = cov(A), weights = smith_hazel, intensity = 1)


###############################################################################
# Desired gains index
###############################################################################


# Calculate desired gains index and resulting trait weights

DG <- make_index(
  genotype.effects         = A,
  var.mat         = V, #posterior variance covariance matrix of BLUPS (ntrait * ngenotype x ntrait * ngenotype), also accepts a simple covariance matric
  desired.gains           = desired_gain
)

# Inspect results
head(DG$index)  # Desired Gains index for each genotype
DG$weights

###############################################################################
# Back-solve average allele-substitution effects
###############################################################################

# Convert breeding values to marker effects via backsolve_marker_effects()
alpha <- backsolve_marker_effects(
  marker.mat       = W,
  G.mat            = GRM,
  genotype.effects = A
)

head(alpha)

###############################################################################
# Strategy: Desired gains; cross expectation (average GEBV)
###############################################################################

# Build all possible parent crosses
PotCrosses <- make_cross_plan(parents = 1:nrow(Mtrain))
head(PotCrosses)

# Get multi-trait cross expectation
expectations <- calc_midparent_inbred(
  crosses         = PotCrosses,
  marker.mat      = Mtrain,
  marker.effects         = alpha,
  weights         = DG$weights,
  nthreads       = 7
)

head(expectations$cross.df)


###############################################################################
# Optimal cross selection (trade-off: diversity vs gain; here using SPV of index)
# Constrain parental contributions
###############################################################################

# Set the maximum number of crosses in which each parent can be used.
# Here, every parent can contribute to at most four crosses.
maximum_contrib <- rep(4, nrow(Mtrain))

# Individual-specific limits can be set by changing the corresponding entry.
# For example, parent 1 can be used at most twice.
maximum_contrib[1] <- 2


# Set the minimum number of crosses in which each parent must be used.
# A value of zero means that the parent is not required to contribute.
minimum_contrib <- rep(0, nrow(Mtrain))

# For example, require parents 2, 3, and 4 to contribute to at least one cross.
minimum_contrib[c(2, 3, 4)] <- 1


# Optimise the crossing plan subject to the parental contribution constraints
ocs_pareto <- optimize_cross_plan(
  candidate.crosses = PotCrosses,
  criterion         = expectations$index.df$GEBV.IDX,
  G.mat             = GRM,
  parents.upper     = maximum_contrib,
  parents.lower     = minimum_contrib,
  method            = "pareto",
  ncrosses          = 50,
  plot              = TRUE
)

# Inspect resulting Pareto solutions
ocs_pareto$pareto.plans[[1]]
head(ocs_pareto$pareto.frontier)
#validate contributions
summary(as.factor(unlist(ocs_pareto$pareto.plans[[1]])))
