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

desired_gains      <- c(10, 10)
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
# Create parents and doubled haploid training population
###############################################################################

parent_pop <- newPop(founder_pop)

rm(founder_pop)

f1_pop <- randCross(
  parent_pop,
  nCrosses = n_crosses,
  simParam = SP
)

dh_pop <- makeDH(
  f1_pop,
  nDH      = n_dh,
  simParam = SP
)

dh_pop <- setPheno(
  dh_pop,
  H2       = h2,
  simParam = SP
)


###############################################################################
# Prepare genomic and phenotypic data
###############################################################################

marker_mat <- pullSnpGeno(dh_pop)

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

p <- colMeans(marker_mat) / 2

# Additive coding: W = M - 2p

W <- sweep(
  marker_mat,
  2,
  2 * p,
  FUN = "-"
)

GRM <- tcrossprod(W) / ncol(marker_mat)


###############################################################################
# Fit multi-trait GBLUP
###############################################################################

asreml.options(
  pworkspace = "5gb",
  workspace  = "5gb"
)

GBLUP_asr <- asreml(
  fixed     = y ~ 1 + trait,
  random    = ~ vm(id, GRM):us(trait),
  residual  = ~ dsum(~ id(units) | trait),
  na.action = na.method(y = "include"),
  data      = pheno_df
)


###############################################################################
# Prediction covariance matrix
###############################################################################

pred <- predict(
  GBLUP_asr,
  classify = "vm(id, GRM):trait",
  only     = "vm(id, GRM):trait",
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

GRM_G <- GRM %x% G

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
rownames(bv_df) <- rownames(marker_mat)


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
# Backsolve marker effects
###############################################################################

alpha <- backsolve_marker_effects(
  marker.mat       = W,
  G.mat            = GRM,
  genotype.effects = bv_df
)

head(alpha)


###############################################################################
# Genetic map
###############################################################################

map_df <- getSnpMap()
map_df$site <- seq_len(nrow(map_df))


###############################################################################
# Candidate crosses
###############################################################################

potential_crosses <- make_cross_plan(
  parents = 1:nrow(marker_mat)
)

head(potential_crosses)


###############################################################################
# Cross expectation
###############################################################################

expectations <- calc_midparent_inbred(
  crosses        = potential_crosses,
  marker.mat     = marker_mat,
  marker.effects = alpha,
  weights        = dg_index$weights,
  p              = p,
  nthreads       = n_threads
)

head(expectations$cross.df)


###############################################################################
# Optimal haploid value
###############################################################################

ohv <- calc_optimal_haploid_value(
  crosses        = potential_crosses,
  marker.mat     = marker_mat,
  marker.effects = alpha,
  weights        = dg_index$weights,
  nthreads       = n_threads
)

head(ohv$cross.df)


###############################################################################
# Segregation variance and superior progeny value
###############################################################################

selection_proportion <- n_crosses / nrow(potential_crosses)

intensity <- dnorm(
  qnorm(1 - selection_proportion)
) / selection_proportion

spv <- calc_spv_inbred(
  crosses        = potential_crosses,
  genetic.map    = map_df,
  marker.mat     = marker_mat,
  marker.effects = alpha,
  t              = 0,
  intensity      = intensity,
  weights        = dg_index$weights,
  p              = p,
  type           = "DH",
  covariance     = FALSE,
  method         = 2,
  nthreads       = n_threads
)

head(spv$cross.df)
head(spv$index.df)


###############################################################################
# Optimal cross selection
###############################################################################

# Pareto optimisation returns alternative gain-diversity trade-offs.

ocs_pareto <- optimize_cross_plan(
  candidate.crosses = potential_crosses,
  criterion         = spv$index.df$SPV.IDX,
  G.mat             = GRM,
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
  candidate.crosses = potential_crosses,
  criterion         = spv$index.df$SPV.IDX,
  G.mat             = GRM,
  method            = "angle",
  ncrosses          = n_crosses,
  target.angle      = 15
)

head(ocs_angle)


###############################################################################
# Select crossing plan
###############################################################################

# Select one Pareto solution according to the breeding objective.

selected_plan_id <- 4000

cross_plan <- ocs_pareto$pareto.plans[[selected_plan_id]]


###############################################################################
# Evaluate crossing plan
###############################################################################

ohv_cols <- paste0("OHV.", trait_names)

cross_metrics <- cbind(
  spv$cross.df,
  ohv$cross.df[ohv_cols]
)

summarize_cross_plan(
  cross.plan = cross_plan,
  cross.df   = cross_metrics
)

plot_cross_plan(
  cross.plan = cross_plan,
  cross.df   = cross_metrics
)
