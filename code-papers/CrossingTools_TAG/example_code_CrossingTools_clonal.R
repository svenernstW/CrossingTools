###############################################################################
# Libraries
###############################################################################

library(AlphaSimR)
library(asreml)
library(CrossingTools)


###############################################################################
# Simulation settings
###############################################################################
n_chr               <- 10
n_parent            <- 50
n_snp_per_chr       <- 500
n_qtl_per_chr       <- 500
n_seg_site_per_chr <- n_snp_per_chr + n_qtl_per_chr
n_crosses           <- 50
n_progeny_per_cross <- 20

# Trait
# Loosely following Werner et al. (2023)
# https://doi.org/10.1007/s00122-023-04300-6

trait_mean <- 0
trait_var  <- 1
mean_dd    <- 0.3
var_dd     <- 0.2
h2         <- 0.5



# Computation
n_threads <- 7




###############################################################################
# Create founders
###############################################################################


founder_pop <- runMacs(
  nInd     = n_parent,
  nChr     = n_chr,
  segSites = n_seg_site_per_chr,
  inbred   = FALSE,
  ploidy   = 2
)


###############################################################################
# Simulation parameters
###############################################################################

SP <- SimParam$new(founder_pop)

SP$addSnpChip(
  nSnpPerChr = n_snp_per_chr
)


SP$addTraitAD(
  nQtlPerChr = n_qtl_per_chr,
  mean       = trait_mean,
  var        = trait_var,
  meanDD     = mean_dd,
  varDD      = var_dd,
  useVarA    = TRUE
)


###############################################################################
# Create parents and training population
###############################################################################

parent_pop <- newPop(founder_pop)

rm(founder_pop)

f1_pop <- randCross(
  parent_pop,
  nCrosses = n_crosses,
  nProgeny = n_progeny_per_cross,
  simParam = SP
)

training_pop <- setPheno(
  f1_pop,
  H2       = h2,
  simParam = SP
)


###############################################################################
# Prepare genomic and phenotypic data
###############################################################################

marker_mat <- pullSnpGeno(training_pop)

# Phased haplotypes are required for segregation variances and SPV.
# In practice, genotype data could be phased using software such as BEAGLE.

haplo_mat <- pullSnpHaplo(training_pop)

hap_1 <- haplo_mat[
  seq(1, 2 * nrow(marker_mat), by = 2),
  ,
  drop = FALSE
]

hap_2 <- haplo_mat[
  seq(2, 2 * nrow(marker_mat), by = 2),
  ,
  drop = FALSE
]

rownames(hap_1) <- rownames(marker_mat)
rownames(hap_2) <- rownames(marker_mat)

pheno_df <- data.frame(
  id = factor(training_pop@id),
  y  = training_pop@pheno[, 1]
)


###############################################################################
# Additive genomic relationship matrix
###############################################################################

p <- colMeans(marker_mat) / 2
q <- 1 - p

# Additive coding: W = M - 2p
W <- sweep(
  marker_mat,
  2,
  2 * p,
  FUN = "-"
)

GRM_A <- tcrossprod(W) / ncol(marker_mat)
GRM_A <- GRM_A + diag(1e-8, nrow(GRM_A))


###############################################################################
# Dominance genomic relationship matrix
###############################################################################

# Heterozygosity indicator
H <- 1 * (marker_mat == 1)

# Statistical dominance-deviation coding:
# Z = H - (1 - 2p)W - 2pq

Z <- H - sweep(
  W,
  2,
  1 - 2 * p,
  FUN = "*"
)

Z <- sweep(
  Z,
  2,
  2 * p * q,
  FUN = "-"
)

GRM_D <- tcrossprod(Z) / ncol(marker_mat)
GRM_D <- GRM_D + diag(1e-8, nrow(GRM_D))


###############################################################################
# Fit additive-dominance GBLUP
###############################################################################

asreml.options(
  pworkspace = "5gb",
  workspace  = "5gb"
)

GBLUP_asr <- asreml(
  fixed     = y ~ 1,
  random    = ~ vm(id, GRM_A) + vm(id, GRM_D),
  residual  = ~ id(units),
  na.action = na.method(y = "include"),
  data      = pheno_df
)


###############################################################################
# Extract breeding values and dominance deviations
###############################################################################

random_coef_df <- as.data.frame(
  GBLUP_asr$coefficients$random
)

bv_idx <- grep(
  "GRM_A",
  rownames(random_coef_df),
  fixed = TRUE
)

dd_idx <- grep(
  "GRM_D",
  rownames(random_coef_df),
  fixed = TRUE
)


bv_df <- random_coef_df[
  bv_idx,
  ,
  drop = FALSE
]

dd_df <- random_coef_df[
  dd_idx,
  ,
  drop = FALSE
]

colnames(bv_df) <- "trait1"
colnames(dd_df) <- "trait1"

rownames(bv_df) <- rownames(marker_mat)
rownames(dd_df) <- rownames(marker_mat)


###############################################################################
# Backsolve marker effects
###############################################################################

alpha <- backsolve_marker_effects(
  marker.mat       = W,
  G.mat            = GRM_A,
  genotype.effects = bv_df
)

delta <- backsolve_marker_effects(
  marker.mat       = Z,
  G.mat            = GRM_D,
  genotype.effects = dd_df
)

head(alpha)
head(delta)


###############################################################################
# Genetic map
###############################################################################

map_df <- getSnpMap()
map_df$site <- seq_len(nrow(map_df))


###############################################################################
# Candidate crosses
###############################################################################

potential_crosses <- make_cross_plan(
  parents = 1:nrow(marker_mat))

head(potential_crosses)


###############################################################################
# Cross expectation
###############################################################################

expectations <- calc_midparent_outcross(
  crosses          = potential_crosses,
  marker.mat       = marker_mat,
  marker.effects.A = alpha,
  marker.effects.D = delta,
  p                = p,
  nthreads         = n_threads
)

head(expectations)


###############################################################################
# Optimal haploid value
###############################################################################

ohv <- calc_optimal_haploid_value(
  crosses        = potential_crosses,
  marker.mat     = marker_mat,
  marker.effects = alpha,
  nthreads       = n_threads
)

head(ohv)


###############################################################################
# Segregation variance and superior progeny value
###############################################################################

selection_proportion <- n_crosses / nrow(potential_crosses)

intensity <- dnorm(
  qnorm(1 - selection_proportion)
) / selection_proportion

spv <- calc_spv_outcross(
  crosses          = potential_crosses,
  genetic.map      = map_df,
  hap.mat1         = hap_1,
  hap.mat2         = hap_2,
  marker.effects.A = alpha,
  marker.effects.D = delta,
  intensity        = intensity,
  p                = p,
  covariance       = FALSE,
  nthreads         = n_threads
)

head(spv)


###############################################################################
# Optimal cross selection
###############################################################################

# Pareto optimisation returns alternative gain-diversity trade-offs.

ocs_pareto <- optimize_cross_plan(
  candidate.crosses = potential_crosses,
  criterion         = spv$TSPV.trait1,
  G.mat             = GRM_A,
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
  criterion         = spv$TSPV.trait1,
  G.mat             = GRM_A,
  method            = "angle",
  ncrosses          = n_crosses,
  target.angle      = 15
)

head(ocs_angle)


###############################################################################
# Select crossing plan
###############################################################################

# Select one Pareto solution according to the breeding objective.

selected_plan_id <- 3000

cross_plan <- ocs_pareto$pareto.plans[[selected_plan_id]]


###############################################################################
# Evaluate crossing plan
###############################################################################

cross_metrics <- cbind(
  spv,
  ohv["OHV.trait1"]
)

summarize_cross_plan(
  cross.plan = cross_plan,
  cross.df   = cross_metrics
)

plot_cross_plan(
  cross.plan = cross_plan,
  cross.df   = cross_metrics
)
