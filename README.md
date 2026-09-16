# CrossingTools

**CrossingTools** is an R package for scalable cross evaluation and mating-plan optimisation in plant breeding. It integrates progeny prediction, multi-trait selection, and genetic-diversity management within a unified framework for supporting crossing decisions.

CrossingTools is designed for applications in line, hybrid, and clonal breeding programmes and uses marker genotypes together with estimated genetic effects to evaluate prospective crosses and optimise crossing plans.

## Main features

CrossingTools provides tools to:

* calculate expected progeny performance for crosses;
* predict within-family segregation variances and covariances for different parental and progeny types;
* evaluate crosses using the expected cross mean, Superior Progeny Value (SPV), and Optimal Haploid Value (OHV);
* account for additive and, where applicable, dominance effects;
* define multi-trait breeding objectives using **Smith–Hazel** economic-weight indices or **Desired Gains** indices;
* optimise mating plans while balancing predicted cross performance and genomic diversity;
* perform optimal cross selection using either a target-angle approach or Pareto-based optimisation; and
* incorporate practical constraints such as mandatory or excluded crosses and restrictions on parental contributions.

Computationally intensive components are implemented in C++, allowing large numbers of candidate crosses to be evaluated efficiently.

## Installation

The current development version can be installed directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("svenernstW/CrossingTools")
```

The package can then be loaded with:

```r
library(CrossingTools)
```

## Documentation

A package manual and example scripts are available in this repository. The example scripts demonstrate the main workflows for cross evaluation, multi-trait index construction, and mating-plan optimisation.

## Reporting problems and requesting features

If you encounter a bug, unexpected behaviour, or another technical problem, please open an issue in the GitHub **Issues** section of this repository. When possible, include a minimal reproducible example together with the relevant input structure, error message, and your R and CrossingTools versions.

Feature requests and suggestions for improving the package are also welcome through GitHub Issues. Keeping these requests on GitHub makes them visible and easier to track.

For collaboration enquiries, questions that are not suitable for a public issue, or other requests, please contact:

**Sven Ernst Weber**
Department of Plant Breeding
Justus Liebig University Giessen
Email: [Sven.Weber@agrar.uni-giessen.de](mailto:Sven.Weber@agrar.uni-giessen.de)

## Acknowledgements

I would particularly like to thank Matthias Frisch for his general guidance and insights on plant breeding, selection theory, and methodology. I also gratefully acknowledge Carola Zenke-Philippi, Eva Herzog, Philipp Heilmann, and Joshua Okoye for many fruitful discussions that contributed to the development of CrossingTools. I further thank Brian Kinghorn and Andrew Kinghorn at the University of New England for valuable discussions on gain–diversity trade-offs and their applications in animal breeding.

## Reference

The theoretical framework and implementation of CrossingTools are described in:

> Weber, S.E., Waters, D.L., Werner, C.R. & Tolhurst, D.J. *CrossingTools: Scalable multi-trait cross evaluation and mating-plan optimisation in plant breeding.* In preparation.

The full citation will be updated following publication.

## License

CrossingTools is distributed under the **GNU General Public License version 3 (GPL-3.0)**. See the `LICENSE` file for details.
