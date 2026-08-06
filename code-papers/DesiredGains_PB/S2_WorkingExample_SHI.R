
#######################################
#
# Supplementary Script S2
#
# Working example:
# Constructing a Smith-Hazel Index
#
# Script author: D.J. Tolhurst
#
#######################################

# This working example constructs a Smith-Hazel Index for
# two genetically correlated traits.
#
# The example uses the simulated data from the previous script,
# where phenotypes were simulated for 2000 inbred genotypes, with
# a genomic relationship matrix constructed from 6000
# markers. The trait heritabilities were 0.7 and 0.4, while
# the genetic correlation was 0.4.
#
# The example involves five steps for constructing the
# Smith-Hazel Index:
#
#  1. Define the breeding objective, as represented by the
#     economic weights vector, a
#  2. Fit a multivariate mixed model to the phenotypic data,
#     producing estimates of the genetic and residual
#     variance-covariance matrices
#  3. Predict multivariate EBVs and compute their
#     variance-covariance matrix
#  4. Construct the Smith-Hazel Index for each selection candidate
#  5. Assess the reliability of the index and select the top
#     candidates as parents for the next generation
#
# The desired gains implied by the economic weights are
# also calculated.

library(AlphaSimR)
library(asreml)


#######################################
# Step 1: Define the breeding objective

# Vector of economic weights
(a <- c(T1 = 1, T2 = 0))
# Assign economic value to improvement in trait 1, with no
# direct economic value assigned to trait 2

# Vector of true aggregate genotypes implied by the economic weights vector
H_SH <- bvs %*% a
names(H_SH) <- rownames(H_SH)


#######################################
# Step 2: Fit a multivariate mixed model to the phenotypic data

str(pheno_df)
asr_shi <- asreml(y ~ trait,
                  random = ~ us(trait):vm(genotype, Gg),
                  residual = ~ dsum(~genotype|trait), # assuming independence between residuals based on the simulation, alter as required
                  data = pheno_df,
                  workspace = 6e8)
asr_shi <- update(asr_shi)
asr_shi <- update(asr_shi)
summary(asr_shi)$varcom


#######################################
# Step 3: Predict multivariate EBVs and compute their
# variance-covariance matrix

pred_shi <- predict(asr_shi,
                    classify = "trait:genotype",
                    only = "us(trait):vm(genotype, Gg)", # depending on asreml version, consider replacing "only" with "onlyuse"
                    vcov = TRUE,
                    pworkspace = 8e8, maxit = 1)

# a) Multivariate genomic estimated breeding values (GEBVs)
gebvs <- matrix(pred_shi$pvals$predicted.value, ncol = ntraits)
colnames(gebvs) <- colnames(G)
rownames(gebvs) <- geno_names
range(c(gebvs) - asr_shi$coefficients$random) # check
plot(gebvs)


# b) Variance-covariance matrix of the multivariate GEBVs
pev_gebvs <- as.matrix(pred_shi$vcov)
range(diag(pev_gebvs) - asr_shi$vcoeff$random) # check
var_gebvs <- kronecker(G, Gg) - pev_gebvs


# c) Calculate the marginal variance-covariance matrix of the GEBVs
# (used for reliability calculations)
V_bar <- matrix(0, ntraits, ntraits)
ii <- split(seq_len(ngenos * ntraits),
            rep(seq_len(ngenos), times = ntraits))
V_ii <- lapply(ii, function(j) var_gebvs[j, j, drop = FALSE])
V_bar <- Reduce("+", V_ii) / ngenos
dimnames(V_bar) <- dimnames(G)
V_bar


# d) Calculate the vector of index coefficients
b_SH <- a
names(b_SH) <- names(d)
b_SH


#######################################
# Step 4: Construct the Smith-Hazel Index for each selection candidate

I_SH <- gebvs %*% b_SH
rownames(I_SH) <- names(I_SH) <- geno_names
hist(I_SH)
plot(I_SH, I_DG); abline(a=0, b=1) # compare to Desired Gains Index


#######################################
# Step 5: Assess the reliability of the index and select the
# top candidates

# a) Obtain model-based and empirical reliability

# REML estimate of the genetic variance-covariance matrix
# obtain REML estimate of the genetic variance-covariance matrix
G_reml <- diag(ntraits)
G_reml[upper.tri(G_reml, diag = TRUE)] <- asr_dgi$vparameters[grep("genotype", names(asr_dgi$vparameters))]
G_reml[2,1] <- G_reml[1,2]
dimnames(G_reml) <- dimnames(G)
G_reml

# Model-based reliability
var_I <- t(b_SH) %*% V_bar %*% b_SH
var_H <- t(a) %*% G_reml %*% a
(rel_model <- var_I / var_H)
# Note that we can also obtain the reliabilities for each selection candidate as
rel_candidate <- sapply(seq_len(ngenos), function(i) t(a) %*% V_ii[[i]] %*% a / (t(a) %*% G_reml %*% a * diag(Gg)[i]))
names(rel_candidate) <- geno_names
hist(rel_candidate)
summary(rel_candidate)

# Empirical reliability, available because the true breeding
# values are known in the simulation
(rel_empirical <- cor(I_SH, H_SH)^2)


# b) Select the top 30 candidates
(selected <- names(sort(I_SH, decreasing = TRUE)[seq_len(nselect)]))
V_bar %*% a * intensity/c(sqrt(var_I)) # expected response in the traits
colMeans(bvs[selected,]) # realised response, available because this is a simulation
plot(gebvs); points(gebvs[selected,], col = "red")
plot(bvs); points(bvs[selected,], col = "red")
mean(H_SH[selected]) # response in the aggregate genotype

# compare to using the true aggregate genotype
(selected_true <- names(sort(H_SH, decreasing = TRUE)[seq_len(nselect)]))
colMeans(bvs[selected_true,])
plot(bvs); points(bvs[selected_true,], col = "red")


#######################################
# Obtaining the implied desired gains

# Desired genetic gains which rank individuals equivalent to
# the Smith-Hazel Index can be obtained by equating the
# Desired Gains Index coefficients to the Smith-Hazel Index
# coefficients. With GEBVs, the vector of index coefficients
# for the Desired Gains Index is equal to the inverse of the
# marginal variance-covariance matrix multiplied by the desired
# gains vector, d

d_implied <- V_bar %*% a
names(d_implied) <- names(a)
d_implied

# Express the implied response relative to trait 1
d_implied <- d_implied / d_implied["T1"]
d_implied


#######################################
# end of script

save.image()

# Save data
save(sim_pop, pheno_df, asr_shi, gebvs, a, V_bar, G, G_reml, b_SH, I_SH, rel_candidate, rel_model, rel_empirical, selected, file = "ScriptS2_SHI.RData")



