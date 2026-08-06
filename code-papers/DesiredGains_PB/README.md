# Supplementary scripts for the Desired Gains and Smith–Hazel indices

This repository contains the supplementary R scripts accompanying the manuscript:

> Werner, C.R., Weber, S.E., Gardner, K.A. & Tolhurst, D.J.  
> *Reviving the Desired Gains Index: An optimal solution for parent selection in public plant breeding programs.*

The scripts provide a working example for constructing a **Desired Gains Index** and a **Smith–Hazel Index** from multivariate genomic estimated breeding values (GEBVs).

## Scripts

### Supplementary Script S1: Desired Gains Index

Script S1 simulates:

- 2,000 inbred genotypes;
- two traits with heritabilities of 0.7 and 0.4;
- a genetic correlation of 0.4;
- 6,000 QTN and 6,000 SNP markers.

It then:

1. defines the desired gains vector;
2. fits a multivariate genomic mixed model;
3. predicts multivariate GEBVs;
4. calculates the marginal variance–covariance matrix of the GEBVs;
5. constructs the Desired Gains Index;
6. calculates index reliability;
7. selects the best 30 candidates;
8. derives equivalent economic weights.

The script saves its main objects to:

```text
ScriptS1_DGI.RData
```

### Supplementary Script S2: Smith–Hazel Index

Script S2 uses the simulated data and objects created in Script S1. It:

1. defines economic weights;
2. fits the multivariate genomic mixed model;
3. predicts multivariate GEBVs;
4. constructs the Smith–Hazel Index;
5. calculates index reliability;
6. selects the best 30 candidates;
7. derives the desired gains implied by the economic weights.

The script saves its main objects to:

```text
ScriptS2_SHI.RData
```

## Requirements

The scripts require R and the following packages:

```r
library(AlphaSimR)
library(asreml)
```
