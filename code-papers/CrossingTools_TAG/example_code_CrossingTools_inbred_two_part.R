###############################################################################
# Libraries
###############################################################################

library(AlphaSimR)
library(asreml)
library(CrossingTools)


###############################################################################
# Simulation settings
###############################################################################

n_chr              <- 10
n_parent           <- 50
n_snp_per_chr      <- 500
n_qtl_per_chr      <- 500
n_seg_site_per_chr <- n_snp_per_chr + n_qtl_per_chr
n_training_crosses <- 50
n_dh               <- 40
n_crosses          <- 50


# Traits
trait_names <- c("trait1", "trait2")
trait_mean  <- c(0, 0)
trait_var   <- c(1, 2)

cor_a <- matrix(
  c(
    1.0, -0.5,
    -0.5,  1.0
  ),
  nrow = 2
)

h2 <- c(0.5, 0.5)

# Selection index

desired_gains       <- c(10, 10)
smith_hazel_weights <- c(1, 1)


# Computation

n_threads <- 7


###############################################################################
# Create founders
###############################################################################

founder_pop <- runMacs(
  nInd     = n_parent,
  nChr     = n_chr,
  segSites = n_seg_site_per_chr,
  inbred   = TRUE,
  species  = "wheat"
)


###############################################################################
# Simulation parameters
###############################################################################

SP <- SimParam$new(founder_pop)

SP$addSnpChip(
  nSnpPerChr = n_snp_per_chr
)

SP$addTraitA(
  nQtlPerChr = n_qtl_per_chr,
  mean       = trait_mean,
  var        = trait_var,
  corA       = cor_a
)


###############################################################################
# Part 1: doubled haploid training population
###############################################################################

parent_pop <- newPop(founder_pop)

rm(founder_pop)

f1_pop <- randCross(
  parent_pop,
  nCrosses = n_training_crosses,
  simParam = SP
)

dh_pop <- makeDH(
  f1_pop,
  nDH      = n_dh,
  simParam = SP
)

# Phenotype DH lines to train genomic prediction models.

dh_pop <- setPheno(
  dh_pop,
  H2       = h2,
  simParam = SP
)


###############################################################################
# Prepare genomic and phenotypic data
###############################################################################

marker_mat_training <- pullSnpGeno(dh_pop)

pheno_df <- data.frame(
  id = factor(
    rep(dh_pop@id, times = length(trait_names))
  ),
  trait = factor(
    rep(trait_names, each = nInd(dh_pop)),
    levels = trait_names
  ),
  y = as.vector(dh_pop@pheno)
)


###############################################################################
# Additive genomic relationship matrix
###############################################################################

p_training <- colMeans(marker_mat_training) / 2

# Centre marker dosages using allele frequencies from the training set.

W_training <- sweep(
  marker_mat_training,
  2,
  2 * p_training,
  FUN = "-"
)

GRM_training <- tcrossprod(W_training) / ncol(marker_mat_training)

# Add a small ridge penalty to ensure positive definiteness

ridge_training <- 1e-4 * mean(diag(GRM_training))

GRM_training <- GRM_training + diag(ridge_training, nrow(GRM_training))

###############################################################################
# Fit multi-trait GBLUP
###############################################################################

asreml.options(
  pworkspace = "5gb",
  workspace  = "5gb"
)

GBLUP_asr <- asreml(
  fixed     = y ~ 1 + trait,
  random    = ~ vm(id, GRM_training):us(trait),
  residual  = ~ dsum(~ id(units) | trait),
  na.action = na.method(y = "include"),
  data      = pheno_df
)


###############################################################################
# Prediction covariance matrix
###############################################################################

pred <- predict(
  GBLUP_asr,
  classify = "vm(id, GRM_training):trait",
  only     = "vm(id, GRM_training):trait",
  vcov     = TRUE
)

# Prediction-error covariance matrix of the breeding-value BLUPs

PEV <- as.matrix(pred$vcov)

# Genetic covariance matrix between traits

G <- diag(
  GBLUP_asr$vparameters[c(1, 3)],
  nrow = 2
)

G[1, 2] <- G[2, 1] <- GBLUP_asr$vparameters[2]

# Prior covariance matrix of the stacked breeding values

GRM_G <- GRM_training %x% G

# Covariance matrix of predicted breeding values:
# Var(A_hat) = Var(A) - PEV

V <- GRM_G - PEV


###############################################################################
# Extract breeding values
###############################################################################

random_coef_df <- as.data.frame(
  GBLUP_asr$coefficients$random
)

bv_df <- cbind(
  random_coef_df[grep(trait_names[1], rownames(random_coef_df), fixed = TRUE), 1],
  random_coef_df[grep(trait_names[2], rownames(random_coef_df), fixed = TRUE), 1]
)

bv_df <- as.data.frame(bv_df)

colnames(bv_df) <- trait_names
rownames(bv_df) <- rownames(marker_mat_training)


###############################################################################
# Evaluate multi-trait selection indices
###############################################################################

# Desired gains specify a 1:1 direction of improvement.

predict_response(
  var.mat       = V,
  desired.gains = desired_gains,
  intensity     = 1
)

# Fixed Smith-Hazel index coefficients.

predict_response(
  var.mat   = V,
  weights   = smith_hazel_weights,
  intensity = 1
)


###############################################################################
# Desired gains index
###############################################################################

dg_index <- make_index(
  genotype.effects = bv_df,
  var.mat          = V,
  desired.gains    = desired_gains
)

head(dg_index$index)
dg_index$weights


###############################################################################
# Smith-Hazel index
###############################################################################

sh_index <- make_index(
  genotype.effects = bv_df,
  var.mat          = V,
  weights          = smith_hazel_weights
)

head(sh_index$index)
sh_index$weights


###############################################################################
# Backsolve marker effects
###############################################################################

alpha <- backsolve_marker_effects(
  marker.mat       = W_training,
  G.mat            = GRM_training,
  genotype.effects = bv_df
)

head(alpha)


###############################################################################
# Part 2: rapid-cycling population
###############################################################################

# Two random mating rounds as an example for the population-improvement part.
# Unlike the DH training population, these candidates are not phenotyped.

rc1_pop <- randCross(
  dh_pop,
  nCrosses = n_crosses,
  nProgeny = 1,
  simParam = SP
)

rc2_pop <- randCross(
  rc1_pop,
  nCrosses = n_crosses,
  nProgeny = 10,
  simParam = SP
)


# RC2 parents may be heterozygous: dosages lose haplotype phase.
# Phased haplotypes preserve linkage needed for segregation variances.

marker_mat_rc2 <- pullSnpGeno(rc2_pop)
marker_hap_rc2 <- pullSnpHaplo(rc2_pop)

# Use training allele frequencies for consistent genomic relationships.

W_rc2 <- sweep(
  marker_mat_rc2,
  2,
  2 * p_training,
  FUN = "-"
)

GRM_rc2 <- tcrossprod(W_rc2) / ncol(W_rc2)

# Add a small ridge penalty to ensure positive definiteness

ridge_rc2 <- 1e-4 * mean(diag(GRM_rc2))

GRM_rc2 <- GRM_rc2 + diag(ridge_rc2, nrow(GRM_rc2))


###############################################################################
# Genetic map
###############################################################################

map_df <- getSnpMap()
map_df$site <- seq_len(nrow(map_df))


###############################################################################
# Candidate crosses
###############################################################################

candidate_crosses <- make_cross_plan(
  parents = 1:nrow(marker_mat_rc2)
)

head(candidate_crosses)

# Each RC2 individual has two rows in marker_hap_rc2 (2i - 1, 2i).
# Use four indices to specify both parental haplotype pairs for SPV.

candidate_crosses_hap <- data.frame(
  p1_h1 = 2 * candidate_crosses[, 1] - 1,
  p1_h2 = 2 * candidate_crosses[, 1],
  p2_h1 = 2 * candidate_crosses[, 2] - 1,
  p2_h2 = 2 * candidate_crosses[, 2]
)

head(candidate_crosses_hap)


###############################################################################
# Cross expectation
###############################################################################

# Apply marker effects and index weights estimated in the DH training set.

expectations <- calc_midparent_inbred(
  crosses        = candidate_crosses,
  marker.mat     = marker_mat_rc2,
  marker.effects = alpha,
  weights        = sh_index$weights,
  p              = p_training,
  nthreads       = n_threads
)

head(expectations$cross.df)
head(expectations$index.df)


###############################################################################
# Optimal haploid value
###############################################################################

ohv <- calc_optimal_haploid_value(
  crosses        = candidate_crosses,
  marker.mat     = marker_mat_rc2,
  marker.effects = alpha,
  weights        = sh_index$weights,
  nthreads       = n_threads
)

head(ohv$cross.df)
head(ohv$index.df)


###############################################################################
# Segregation variance and superior progeny value
###############################################################################

selection_proportion <- n_crosses / nrow(candidate_crosses)

intensity <- dnorm(
  qnorm(1 - selection_proportion)
) / selection_proportion

spv <- calc_spv_inbred(
  crosses        = candidate_crosses_hap,
  genetic.map    = map_df,
  marker.mat     = marker_hap_rc2,
  marker.effects = alpha,
  t              = 1,
  intensity      = intensity,
  weights        = sh_index$weights,
  p              = p_training,
  type           = "DH",
  covariance     = FALSE,
  method         = 1,
  nthreads       = n_threads
)

head(spv$cross.df)
head(spv$index.df)


###############################################################################
# Optimal cross selection
###############################################################################

# Use RC2 relationships to balance gain and diversity among candidates.
# Pareto optimisation returns alternative gain-diversity trade-offs.

ocs_pareto <- optimize_cross_plan(
  candidate.crosses = candidate_crosses,
  criterion         = spv$index.df$SPV.IDX,
  G.mat             = GRM_rc2,
  method            = "pareto",
  ncrosses          = n_crosses,
  plot              = TRUE
)

ocs_pareto$pareto.plans[[1]]

head(
  ocs_pareto$pareto.frontier
)


# Alternatively, optimise a single gain-diversity trade-off.

ocs_angle <- optimize_cross_plan(
  candidate.crosses = candidate_crosses,
  criterion         = spv$index.df$SPV.IDX,
  G.mat             = GRM_rc2,
  method            = "angle",
  ncrosses          = n_crosses,
  target.angle      = 15
)

head(ocs_angle)


###############################################################################
# Select crossing plan
###############################################################################

# Select one Pareto solution according to the breeding objective.

selected_plan_id <- 1000 #example!

selected_cross_plan <- ocs_pareto$pareto.plans[[selected_plan_id]]


###############################################################################
# Evaluate crossing plan
###############################################################################

ohv_cols <- paste0("OHV.", trait_names)

cross_metrics <- cbind(
  candidate_crosses,
  spv$cross.df[, 5:ncol(spv$cross.df), drop = FALSE],
  ohv$cross.df[, ohv_cols, drop = FALSE]
)

summarize_cross_plan(
  cross.plan = selected_cross_plan,
  cross.df   = cross_metrics
)

plot_cross_plan(
  cross.plan = selected_cross_plan,
  cross.df   = cross_metrics
)

