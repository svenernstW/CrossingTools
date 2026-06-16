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
# Fill breeding pipeline
###############################################################################

# Stage 2: create F1s and doubled haploids (DH)
F1 = randCross(Parents, 50)
DH = makeDH(F1, 40)

# Phenotype DH lines
DH = setPheno(DH,H2 = c(0.5,0.5))

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
# Genomic relationship matrix (GRM)
###############################################################################

# Compute GRM with a simple alternative, doesn't really matter
GRM <- tcrossprod(Mtrain)/ncol(Mtrain)


###############################################################################
# Fit GBLUP model (ASReml) and predict breeding values
###############################################################################

# Increase ASReml workspace 
asreml.options(pworkspace = 1e+9,workspace = 1e+9)

GBLUP <- asreml(fixed = value ~ 1+trait,
                random = ~us(trait):vm(id,GRM),
                residual = ~ dsum(~ id(units) | trait),
                na.action = na.method(y = "include"),
                data = Pheno)




# EBVs (BLUPs) for id x trait, plus vcov matrix of those predictions
pred <- predict(
  GBLUP,
  classify = "id:trait",   # or "id:trait", just be consistent
  vcov = TRUE
)


PEV <- as.matrix(pred$vcov)       # this is your PEV(A)


G <- diag(GBLUP$vparameters[c(1,3)],2,2) #this is the genetic covariance matrix
G[2,1] <- G[1,2] <- GBLUP$vparameters[c(2)]


GRM.G <- GRM %x% G # this is your var(A)

V <- GRM.G-PEV # this is your var(A_tilde)


###############################################################################
# Extract breeding values (two traits)
###############################################################################

# Breeding values for Trait1 and Trait2 (stacked random coefficients)
A <- as.data.frame(cbind(
  GBLUP$coefficients$random[1:nInd(DH)],
  GBLUP$coefficients$random[(nInd(DH)+1):(2*nInd(DH))]
))

names(A) <- c("Trait1","Trait2")

###############################################################################
# Summarize expected index behavior from trait covariance (quick intuition)
###############################################################################
desired_gain <- c(10, 10)  # desired gains per trait
smith_hazel  <- c(1, 1)    # economic weights per trait

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
# Marker effects for the two traits
###############################################################################

# Convert breeding values to marker effects via get_marker_effects()
U_12 <- backsolve_marker_effects(
  marker.mat     = Mtrain,
  G.mat          = GRM,
  genotype.effects        = A
)

head(U_12)
###############################################################################
# Genetic map for markers
###############################################################################

# Get marker genetic map and assign sequential site IDs
map        <- getSnpMap()
map$site   <- seq_len(nrow(map))
#map$site   <- map$id

###############################################################################
# Strategy: Desired gains; cross expectation (average GEBV), Optimal haploid value and superior progeny index
###############################################################################

# Build all possible parent crosses
PotCrosses <- make_cross_plan(parents = 1:nrow(Mtrain))
head(PotCrosses)

# Get multi-trait cross expectation
expectations <- calc_midparent_inbred(
  crosses         = PotCrosses,
  marker.mat      = Mtrain,
  marker.effects         = U_12,
  weights         = DG$weights,
  nthreads       = 7
)

head(expectations$cross.df)

# Get multi-trait cross optimal haploid value
ohv <- calc_optimal_haploid_value(
  crosses         = PotCrosses,
  marker.mat      = Mtrain,
  marker.effects         = U_12,
  weights         = DG$weights,
  nthreads       = 7
)

head(ohv$cross.df)

# Get multi-trait segregation var & covariance and superior progeny value
# Selection intensity at cross level
nCrosses <- 50
alpha     <- nCrosses / nrow(PotCrosses) # this is really just a weighting factor on the segregation variance
intensity <- dnorm(qnorm(1 - alpha), 0, 1) / alpha


#with two traits and about two million crosses as in this example might take a while
spv <- calc_spv_inbred(
  crosses         = PotCrosses,
  genetic.map     = map,
  marker.mat      = Mtrain,
  marker.effects         = U_12,
  t               = 0, # how many rounds of random mating before DH or RIL creation
  intensity       = intensity,
  weights         = DG$weights,
  type = "DH",
  covariance      = F,
  method          = 2, # 1 for lehermeier and 2 for osthushenrich
  nthreads       = 7)

###############################################################################
# Optimal cross selection (trade-off: diversity vs gain; here using SPV of index)
###############################################################################

# Notes:
# - Can add fixed.crosses (always conducted)
# - Can remove potential crosses from optimization by removing them from `crosses`
ocs_pareto <- optimize_cross_plan(
  candidate.crosses       = PotCrosses[-1,],
  fixed.crosses = PotCrosses[1,,drop=F],
  criterion             = spv$index$SPV.IDX[-1],
  criterion.fixed = spv$index$SPV.IDX[1],
  G.mat             = GRM,
  method = "pareto", #either pareto to return a pareto with multiple solutions or "angle" to only find a single solution to maximize gain and balance diversity alon a given target angle
  ncrosses      = 50,
  plot=T
)

ocs_pareto$pareto.plans[[1]] # a list of all pareto plans
head(ocs_pareto$pareto.frontier) # coordinates on the pareto

#alternativly, optimize a long a single target angle between gain and genetic similarity

ocs_angle <- optimize_cross_plan(
  candidate.crosses       = PotCrosses[-1,],
  fixed.crosses = PotCrosses[1,,drop=F],
  criterion             = spv$index$SPV.IDX[-1],
  criterion.fixed = spv$index$SPV.IDX[1],
  G.mat             = GRM,
  method = "angle", #either pareto to return a pareto with multiple solutions or "angle" to only find a single solution to maximize gain and balance diversity alon a given target angle
  ncrosses      = 50,
  target.angle  = 15
)
head(ocs_angle)

###############################################################################
# Final cross plan and evaluation
###############################################################################

# Final cross plan
crosses <- ocs_pareto$pareto.plans[[4950]] # pick the one that fits your goals best


###############################################################################
# Evaluate and plot a crossing plan
###############################################################################
# Evaluate the plan
cross.df <- spv$cross.df

cross.df <- cbind(spv$cross.df,ohv$cross.df[,3:4])


summarize_cross_plan(
  cross.plan  = crosses,
  cross.df = cross.df
)

# Plot the plan
plot_cross_plan(
  cross.plan  = crosses,
  cross.df = cross.df
)

