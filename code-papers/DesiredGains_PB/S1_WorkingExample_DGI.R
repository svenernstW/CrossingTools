
#######################################
#
# Supplementary Script S1
#
# Working example:
# Constructing a Desired Gains Index
#
# Script author: D.J. Tolhurst
#
#######################################

# This working example constructs a Desired Gains Index for
# two genetically correlated traits.
#
# Phenotypes are simulated for 2000 inbred genotypes, with
# a genomic relationship matrix constructed from 6000
# markers. The trait heritabilities are 0.7 and 0.4, while
# the genetic correlation is 0.4.
#
# The example involves an initial simulation step followed by
# five steps for constructing the Desired Gains Index:
#
#  1. Define a common breeding objective, as represented by the
#     desired gains vector, d
#  2. Fit a multivariate mixed model to the phenotypic data,
#     producing estimates of the genetic and residual
#     variance-covariance matrices
#  3. Predict multivariate EBVs and compute their marginal
#     variance-covariance matrix
#  4. Construct the Desired Gains Index for each selection candidate
#  5. Assess the reliability of the index and select the top
#     candidates as parents for the next generation
#
# The economic weights to achieve the desired gains
# are also calculated.

library(AlphaSimR)
library(asreml)


#######################################
# Step 0: Simulate data for the working example

# a) Simulation parameters
ngenos <- 2000
nQTN <- 500
nSNP <- 500
nchromosomes <- 12
nQTN*nchromosomes # 6000
nSNP*nchromosomes # 6000
nselect <- 30
(intensity <- dnorm(qnorm(1-nselect/ngenos))*ngenos/nselect)
ntraits <- 2


# b) Population variance-covariance matrices
h2  <- c(0.7, 0.4) # heritabilities
mean_g <- c(0,0)   # trait means
var_g <- c(1,1)    # genetic variances
cor_g <- 0.4       # genetic correlation
Cg <- matrix(c(1,cor_g,cor_g,1), ncol = ntraits) # genetic correlation matrix
G <- diag(ntraits)
G[1,2] <- G[2,1] <- cor_g

(var_e <- var_g * (1 - h2) / h2) # residual variances
cor_e <- 0 # residual correlation
R <- diag(var_e) # residual correlation matrix

dimnames(G) <- dimnames(R) <- list(c("T1", "T2"), c("T1", "T2"))
G
R
(P <- G + R)


# c) Simulate a population with 2000 individuals
set.seed(123)
pop <- runMacs(nInd = ngenos,
               nChr = nchromosomes,
               segSites = nQTN + nSNP,
               inbred = TRUE,
               species = "GENERIC")
SP <- SimParam$new(pop)

SP$addTraitA(nQTN,
             mean = mean_g,
             var = var_g,
             corA = Cg)
SP$addSnpChip(nSNP)
sim_pop <- newPop(pop)


# d) Construct the genomic relationship matrix
M <- scale(pullSnpGeno(sim_pop), center = TRUE, scale = FALSE)
Gg <- M %*% t(M)
diag(Gg) <- diag(Gg) + 1e-6 # Add a small ridge to ensure positive definiteness
Gg <- Gg/mean(diag(Gg)) # scale

geno_names <- paste0("G", seq_len(ngenos))
dimnames(Gg) <- list(geno_names, geno_names)
Gg[1:10,1:10]
summary(diag(Gg))
hist(rowMeans(abs(Gg)))


# e) Obtain true breeding values for the two traits
bvs <- sim_pop@gv
plot(bvs)
colnames(bvs) <- colnames(G)
rownames(bvs) <- geno_names


# f) Simulate error and generate phenotypes
error <- scale(matrix(rnorm(ngenos * ntraits), ncol = ntraits)) %*% chol(R)
y <- bvs + error
plot(y)
colnames(y) <- colnames(G)
rownames(y) <- geno_names


# construct phenotypic dataframe
pheno_df <- data.frame(trait = factor(paste0("T", rep(1:2, each = ngenos))),
                       genotype = factor(geno_names, levels = geno_names),
                       y = c(y))
head(pheno_df)

# the simulated effects above will be used to construct the
# Desired Gains Index for the slection candidates.


#######################################
# Step 1: Define a common breeding objective

# Vector of desired genetic gains
(d <- c(T1 = 1, T2 = 0))
# Improve trait 1 while maintaining trait 2 at its current value

# Vector of aggregate genotypes implied by the desired gains vector
# based on the true breeding values
H_DG_true <- bvs %*% solve(G) %*% d
names(H_DG_true) <- rownames(H_DG_true)


#######################################
# Step 2: Fit a multivariate mixed model to the phenotypic data

str(pheno_df)
asr_dgi <- asreml(y ~ trait,
                  random = ~ us(trait):vm(genotype, Gg),
                  residual = ~ dsum(~genotype|trait), # assuming independence between residuals based on the simulation, alter as required
                  data = pheno_df,
                  workspace = 6e8)
asr_dgi <- update(asr_dgi)
asr_dgi <- update(asr_dgi)
summary(asr_dgi)$varcom


#######################################
# Step 3: Predict multivariate EBVs and compute their marginal
# variance-covariance matrix

pred_dgi <- predict(asr_dgi,
                    classify = "trait:genotype",
                    only = "us(trait):vm(genotype, Gg)", # depending on asreml version, consider replacing "only" with "onlyuse"
                    vcov = TRUE,
                    pworkspace = 8e8, maxit = 1)

# a) Multivariate genomic estimated breeding values (GEBVs)
gebvs <- matrix(pred_dgi$pvals$predicted.value, ncol = ntraits)
colnames(gebvs) <- colnames(G)
rownames(gebvs) <- geno_names
range(c(gebvs) - asr_dgi$coefficients$random) # check
plot(gebvs)


# b) Variance-covariance matrix of the multivariate GEBVs
pev_gebvs <- as.matrix(pred_dgi$vcov)
range(diag(pev_gebvs) - asr_dgi$vcoeff$random) # check
var_gebvs <- kronecker(G, Gg) - pev_gebvs


# c) Calculate the marginal variance-covariance matrix of the GEBVs
V_bar <- matrix(0, ntraits, ntraits)
ii <- split(seq_len(ngenos * ntraits),
            rep(seq_len(ngenos), times = ntraits))
V_ii <- lapply(ii, function(j) var_gebvs[j, j, drop = FALSE])
V_bar <- Reduce("+", V_ii) / ngenos
dimnames(V_bar) <- dimnames(G)
V_bar


# d) Calculate the common vector of index coefficients
b_DG <- solve(V_bar) %*% d
names(b_DG) <- names(d)
b_DG


#######################################
# Step 4: Construct the Desired Gains Index for each selection candidate

I_DG <- gebvs %*% b_DG
rownames(I_DG) <- names(I_DG) <- geno_names
hist(I_DG)


#######################################
# Step 5: Assess the reliability of the index and select the
# top candidates

# a) Obtain model-based and empirical reliability

# obtain REML estimate of the genetic variance-covariance matrix
G_reml <- diag(ntraits)
G_reml[upper.tri(G_reml, diag = TRUE)] <- asr_dgi$vparameters[grep("genotype", names(asr_dgi$vparameters))]
G_reml[2,1] <- G_reml[1,2]
dimnames(G_reml) <- dimnames(G)
G_reml

# model-based reliability
cov_HI <- t(d) %*% solve(V_bar) %*% d
var_H <- t(d) %*% solve(V_bar) %*% G_reml %*% solve(V_bar) %*% d
var_I <- t(d) %*% solve(V_bar) %*% d # used below for expected response
(rel_model <- cov_HI / var_H)
# Note that we can also obtain the reliabilities for each selection candidate as
rel_candidate <- sapply(seq_len(ngenos), function(i) t(d) %*% solve(V_bar) %*% V_ii[[i]] %*% solve(V_bar) %*% d / (t(d) %*% solve(V_bar) %*% G_reml %*% solve(V_bar) %*% d * diag(Gg)[i]))
names(rel_candidate) <- geno_names
hist(rel_candidate)
summary(rel_candidate)

# empirical reliability, available because the true breeding
# values are known in the simulation
plot(I_DG, H_DG_true)
cor(I_DG, H_DG_true)^2
# Now compare to the true aggregate genotype based on the
# predicted BLUPs with marginal variance-covariance matrix
H_DG <- bvs %*% solve(V_bar) %*% d
names(H_DG) <- rownames(H_DG)
plot(I_DG, H_DG)
(rel_empirical <- cor(I_DG, H_DG)^2)


# b) Select the top 30 candidates
(selected <- names(sort(I_DG, decreasing = TRUE)[seq_len(nselect)]))
d # desired response ratio across traits
d * intensity/c(sqrt(var_I)) # expected response in the traits
colMeans(bvs[selected,]) # realised response, available because this is a simulation
plot(gebvs); points(gebvs[selected,], col = "red")
plot(bvs); points(bvs[selected,], col = "red")

# compare to using the true aggregate genotype
(selected_true <- names(sort(H_DG, decreasing = TRUE)[seq_len(nselect)]))
colMeans(bvs[selected_true,])
plot(bvs); points(bvs[selected_true,], col = "red")


#######################################
# Obtaining the implied economic weights

# Economic weights which rank individuals equivalent to the
# Desired Gains Index can be obtained by equating the
# Smith-Hazel Index coefficients to the Desired Gains Index
# coefficients. With GEBVs, the vector of index coefficients
# for the Smith-Hazel Index is equal to the economic weights, a.

a_equiv <- b_DG
I_SH <- gebvs %*% a_equiv
rownames(I_SH) <- names(I_SH) <- geno_names
hist(I_SH)

plot(I_SH, I_DG); abline(a=0, b=1)



#######################################
# end of script

save.image()

# Save data for use in Script S2, Constructing a Smith-Hazel Index
save(sim_pop, pheno_df, asr_dgi, gebvs, d, V_bar, G, G_reml, b_DG, I_DG, rel_candidate, rel_model, rel_empirical, selected, file = "ScriptS1_DGI.RData")


