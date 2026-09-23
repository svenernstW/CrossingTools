#' Calculate segregation variance and superior progeny values for outcrosses
#'
#' Calculates expected genomic breeding values, expected total genetic values,
#' additive and dominance segregation variances, and superior progeny values
#' for proposed F1 crosses.
#'
#' Additive superior progeny values are calculated from the breeding-value
#' segregation variance, whereas total superior progeny values are calculated
#' from the total genetic segregation variance, including breeding-value
#' variance, dominance-deviation variance, and their covariance.
#'
#' Segregation (co)variances are calculated analytically from phased parental
#' haplotypes, recombination fractions, and additive and dominance marker
#' effects, following the framework of Bonk et al. (2016) for Mendelian
#' sampling covariability among full-sib progeny.
#'
#' Additive genetic values are expressed as breeding values using centred
#' marker genotypes, \eqn{M - 2p}, and average allele-substitution effects
#' (\eqn{\alpha}). Dominance is represented by statistical dominance
#' deviations. The expected total genetic value is therefore the sum of the
#' expected breeding value and expected dominance deviation.
#'
#' For multiple traits, the corresponding segregation covariance matrices are
#' calculated analogously. These quantities can be combined with a selection
#' intensity to obtain additive superior progeny values (SPV) and total
#' superior progeny values (TSPV).
#'
#' @param crosses Matrix or data frame with two columns specifying the parents
#'   of each proposed cross. Parent identifiers may be row indices of
#'   \code{hap.mat1} and \code{hap.mat2}, or character identifiers matching
#'   their row names.
#' @param genetic.map Data frame describing the genetic map. It must contain:
#'   \describe{
#'     \item{\code{site}}{Marker identifier given as a column index or
#'     as a marker name matching \code{colnames(hap.mat1)}.}
#'     \item{\code{chr}}{Chromosome identifier.}
#'     \item{\code{pos}}{Marker position on the chromosome in Morgan.}
#'   }
#'   All markers in the haplotype matrices must be represented.
#' @param hap.mat1 Numeric haplotype matrix with individuals in rows and markers
#'   in columns, containing the first haplotype of each individual.
#' @param hap.mat2 Numeric haplotype matrix with the same dimensions as
#'   \code{hap.mat1}, containing the second haplotype of each individual.
#' @param marker.effects.A Numeric matrix of average allele-substitution effects
#'   (\eqn{\alpha}) with markers in rows and traits in columns. Effects should be
#'   parameterised relative to the reference population defined by \code{p}.
#'   If \code{p = NULL}, the reference allele frequencies are derived from
#'   \code{hap.mat1} and \code{hap.mat2}. Its number of rows must equal
#'   \code{ncol(hap.mat1)}.
#' @param marker.effects.D Numeric matrix of dominance marker effects
#'   (\eqn{\delta}) with markers in rows and traits in columns, used with the
#'   statistical dominance-deviation parameterisation. It must have the same
#'   dimensions as \code{marker.effects.A}.
#' @param intensity Numeric scalar giving the standardized selection intensity
#'   used to calculate SPV and TSPV. The default is 1.
#' @param weights Optional numeric vector with one weight per trait. When
#'   supplied, weighted index values are calculated for each cross.
#' @param p Optional numeric vector of reference allele frequencies, with one
#'   value per marker. If \code{NULL}, allele frequencies are calculated from
#'   \code{hap.mat1} and \code{hap.mat2}. The supplied frequencies are used to
#'   centre marker genotypes for breeding-value calculations.
#' @param covariance Logical. If \code{TRUE}, also calculate breeding-value,
#'   dominance-deviation, and combined additive--dominance segregation
#'   covariance matrices among traits for each cross.
#' @param nthreads Positive integer. Number of computational threads.
#'
#' @return If neither \code{weights} nor trait covariances are requested, a data
#'   frame containing the parental identifiers and, for each trait,
#'   \code{GEBV.<trait>}, \code{TGV.<trait>}, \code{var.A.<trait>},
#'   \code{SPV.<trait>}, \code{var.D.<trait>},
#'   \code{var.TGV.<trait>}, and \code{TSPV.<trait>}.
#'
#'   If \code{weights} is supplied, \code{index.df} contains
#'   \code{GEBV.IDX}, \code{TGV.IDX}, \code{var.A.IDX},
#'   \code{SPV.IDX}, \code{var.D.IDX}, \code{var.TGV.IDX},
#'   and \code{TSPV.IDX}. The total genetic segregation variance
#'   \code{var.TGV.IDX} includes breeding-value variance,
#'   dominance-deviation variance, and their covariance.
#'
#'   If \code{covariance = TRUE}, a list containing:
#'   \describe{
#'     \item{\code{cross.df}}{Cross-specific trait values and segregation
#'     variances.}
#'     \item{\code{index.df}}{Weighted index values, if \code{weights} is
#'     supplied.}
#'     \item{\code{additive.covariances}}{Breeding-value segregation covariance
#'     matrices for each cross.}
#'     \item{\code{dominance.covariances}}{Dominance-deviation segregation
#'     covariance matrices for each cross.}
#'     \item{\code{additive.dominance.covariances}}{Combined additive--dominance
#'     covariance matrices, \eqn{\Sigma_{AD} + \Sigma_{DA}}, for each cross.}
#'   }
#'
#' @references
#' Bonk, S., Reichelt, M., Teuscher, F., Segelke, D. and Reinsch, N. (2016).
#' Mendelian sampling covariability of marker effects and genetic values.
#' \emph{Genetics Selection Evolution}, 48, 36.
#' \doi{10.1186/s12711-016-0214-0}
#'
#' @export

calc_spv_outcross <- function(crosses, genetic.map, hap.mat1, hap.mat2, marker.effects.A, marker.effects.D,
                                   intensity=1, weights = NULL,p=NULL, covariance = FALSE,
                                   nthreads = 4L) {

  traits <- colnames(marker.effects.A)

  if (is.null(traits)) {
    traits <- paste0("trait", seq_len(ncol(marker.effects.A)))
  }

  if(!ncol(crosses) %in% c(2)){stop("ncol(crosses) needs to be 2 ")}
  crosses_in <- crosses

  hap.mat1 <- as.matrix(hap.mat1); hap.mat2 <- as.matrix(hap.mat2)
  effects.A <- as.matrix(marker.effects.A); effects.D <- as.matrix(marker.effects.D)
  crosses <- as.matrix(crosses)

  n.Threads <- nthreads

  if (is.numeric(crosses) || is.integer(crosses)) {
    if (any(!is.finite(crosses))) stop("`crosses` contains non-finite entries.")
    if (any(crosses < 1 | crosses > nrow(hap.mat1))) {
      stop("Some genotype indices in `crosses` are outside 1..nrow(hap.mat1).")
    }
    crosses2 <- as.data.frame(crosses)
  } else {
    if (is.null(rownames(hap.mat1))) stop("Character `crosses` requires rownames(hap.mat1).")

    idx <- match(as.vector(crosses), rownames(hap.mat1))
    if (anyNA(idx)) stop("Some entries in `crosses` are not in rownames(hap.mat1).")

    idx <- matrix(idx, nrow = nrow(crosses), ncol = ncol(crosses), byrow = FALSE)
    crosses2 <- as.data.frame(idx)
  }

  crosses2 <- as.matrix(crosses2)
  storage.mode(crosses2) <- "integer"
  if (anyNA(crosses2)) stop("Internal error: `crosses2` contains NA after conversion.")


  if (!is.logical(covariance) || length(covariance) != 1L) {
    stop("`covariance` must be a single logical (TRUE/FALSE).")
  }
  if (!is.numeric(intensity) || length(intensity) != 1L) {
    warning("`intensity` must be a single numeric (standardized selection differential), setting it to 1.")
    intensity <- 1
  }

  if (length(n.Threads) != 1L || !is.finite(n.Threads) || n.Threads < 1 || n.Threads != as.integer(n.Threads)) {
    stop("`n.Threads` must be a positive integer.")
  }
  nThreads <- as.integer(n.Threads)

  if (ncol(hap.mat1) <= 0L) stop("hap.mat1 must have markers in columns.")
  if (!identical(dim(hap.mat1), dim(hap.mat2))) stop("hap.mat1 and hap.mat2 must have identical dimensions.")
  if (nrow(effects.A) != ncol(hap.mat1)) {
    stop("effects.A must have nrow(effects.A) == ncol(hap.mat1). Found: nrow(effects.A) = ", nrow(effects.A), ", ncol(hap.mat1) = ", ncol(hap.mat1), ".")
  }
  if (nrow(effects.D) != ncol(hap.mat1)) {
    stop("effects.D must have nrow(effects.D) == ncol(hap.mat1). Found: nrow(effects.D) = ", nrow(effects.D), ", ncol(hap.mat1) = ", ncol(hap.mat1), ".")
  }
  if (ncol(effects.D) != ncol(effects.A)) stop("effects.A and effects.D must have the same number of trait columns.")
  if (ncol(crosses) != 2L) stop("`crosses` must have exactly 2 columns (P1, P2).")

  ###############################################################################
  # Reference allele frequencies
  ###############################################################################

  if (is.null(p)) {

    # Genotype dosage matrix M = hap.mat1 + hap.mat2
    p <- colMeans(hap.mat1 + hap.mat2) / 2

  } else {

    if (!is.numeric(p)) {
      stop("`p` must be a numeric vector of reference allele frequencies.")
    }

    p <- as.numeric(p)

    if (length(p) != ncol(hap.mat1)) {
      stop(
        "`p` must contain one allele frequency per marker: ",
        "length(p) = ", length(p),
        ", ncol(hap.mat1) = ", ncol(hap.mat1), "."
      )
    }
  }

  if (any(!is.finite(p))) {
    stop("`p` must contain only finite values.")
  }

  if (any(p < 0 | p > 1)) {
    stop("All values in `p` must be between 0 and 1.")
  }

  ntraits <- ncol(effects.A)

  if (ntraits == 1L && !is.null(weights)) {
    warning("Single trait: weights are ignored (index equals trait).", call. = FALSE)
    weights <- NULL
  }
  if (ntraits == 1L) {
    covariance <- FALSE
  }

  calculate.index <- !is.null(weights) && covariance
  calculate.simple.index <- !is.null(weights) && !covariance

  if (is.null(weights)) {
    weights <- rep(1, ntraits)
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != ntraits) {
      stop("`weights` must have length equal to ncol(effects.A).", call. = FALSE)
    }
    if (any(!is.finite(weights))) {
      stop("`weights` must contain only finite values.", call. = FALSE)
    }
  }

  if (calculate.simple.index) {
    effects.A <- cbind(effects.A, effects.A %*% weights)
    effects.D <- cbind(effects.D, effects.D %*% weights)
  }

  req_cols <- c("site", "chr", "pos")
  if (!all(req_cols %in% names(genetic.map))){
    stop("`genetic.map` must contain columns: site, chr, pos.")
  }
  site <- genetic.map$site

  # Only needed if site is character
  if (is.numeric(site) || is.integer(site)) {

    site_idx <- as.integer(site)

    if (any(site_idx < 1 | site_idx > ncol(hap.mat1), na.rm = TRUE))
      stop("`genetic.map$site` contains indices outside 1..ncol(hap.mat1).")

    if (!all(seq_len(ncol(hap.mat1)) %in% site_idx))
      stop("Some markers in hap.mat1 and 2 are missing from `genetic.map$site`.")

  } else {

    mnames <- colnames(hap.mat1)
    if (is.null(mnames))
      stop("Character `genetic.map$site` requires colnames(hap.mat1).")

    site_chr <- as.character(site)

    if (anyDuplicated(site_chr))
      stop("`genetic.map$site` contains duplicated marker names.")

    if (!all(mnames %in% site_chr))
      stop("Some markers in hap.mat1 are missing from `genetic.map$site`.")

    site_idx <- match(site_chr, mnames)
    if (anyNA(site_idx))
      stop("`genetic.map$site` contains names not found in colnames(hap.mat1).")
  }

  # Order map by chr/pos and carry indices along
  o    <- order(genetic.map$chr, genetic.map$pos)
  map2 <- genetic.map[o, , drop = FALSE]
  ord  <- site_idx[o]


  # relabel sites to 1..p for downstream code
  map2$site <- seq_len(ncol(hap.mat1))


  hap.mat1  <- hap.mat1[, ord, drop = FALSE]
  hap.mat2  <- hap.mat2[, ord, drop = FALSE]
  effects.A <- effects.A[ord, , drop = FALSE]
  effects.D <- effects.D[ord, , drop = FALSE]
  p         <- p[ord]

  chr_levels <- unique(map2$chr)
  genmap_list <- lapply(chr_levels, function(cc) {
    as.matrix(map2$pos[map2$chr == cc])
  })

  temp <- cpp_calculate_covariance_wolfe(
    Crosses    = crosses2,
    genMap     = genmap_list,
    Hap1       = hap.mat1,
    Hap2       = hap.mat2,
    U          = effects.A,
    D          = effects.D,
    intensity  = intensity,
    weights    = weights,
    p          = p,
    covariance = covariance,
    calcindex  = calculate.index,
    nThreads   = nThreads
  )


  crosses_df <- as.data.frame(crosses_in, stringsAsFactors = FALSE)
  names(crosses_df) <- c("parent1", "parent2")
  name_vec <- c(
    paste0("GEBV.", traits),
    paste0("TGV.", traits),
    paste0("var.A.", traits),
    paste0("SPV.", traits),
    paste0("var.D.", traits),
    paste0("var.TGV.", traits),
    paste0("TSPV.", traits)
  )

  crosses_df <- as.data.frame(crosses_in, stringsAsFactors = FALSE)
  names(crosses_df) <- c("parent1", "parent2")

  if (covariance) {
    cv <- as.data.frame(temp$cross_values)

    temp1 <- cv[, 1:(7 * ntraits), drop = FALSE]
    names(temp1) <- name_vec
    temp1 <- cbind(crosses_df, temp1)

    if (calculate.index) {
      temp2 <- cv[, (7 * ntraits + 1):(7 * ntraits + 7), drop = FALSE]
      names(temp2) <- c(
        "GEBV.IDX",
        "TGV.IDX",
        "var.A.IDX",
        "SPV.IDX",
        "var.D.IDX",
        "var.TGV.IDX",
        "TSPV.IDX"
      )
      temp2 <- cbind(crosses_df, temp2)

      if (isTRUE(temp$check_psd)) {
        warning("Some segregation covariance matrices were not psd; for index calculation the nearest psd was used in these places.", call. = FALSE)
      }

      return(list(
        cross.df = temp1,
        index.df = temp2,
        additive.covariances = temp$covA,
        dominance.covariances = temp$covD,
        additive.dominance.covariances = temp$covAD
      ))
    } else {
      return(list(
        cross.df = temp1,
        additive.covariances = temp$covA,
        dominance.covariances = temp$covD,
        additive.dominance.covariances = temp$covAD
      ))
    }
  }

  if (calculate.simple.index) {
    cv <- as.data.frame(temp)

    ntraits.cpp <- ncol(effects.A)  # original traits + 1 appended index trait

    trait_cols <- c(
      seq_len(ntraits),
      ntraits.cpp + seq_len(ntraits),
      2 * ntraits.cpp + seq_len(ntraits),
      3 * ntraits.cpp + seq_len(ntraits),
      4 * ntraits.cpp + seq_len(ntraits),
      5 * ntraits.cpp + seq_len(ntraits),
      6 * ntraits.cpp + seq_len(ntraits)
    )

    idx_cols <- c(
      ntraits.cpp,
      2 * ntraits.cpp,
      3 * ntraits.cpp,
      4 * ntraits.cpp,
      5 * ntraits.cpp,
      6 * ntraits.cpp,
      7 * ntraits.cpp
    )

    temp1 <- cv[, trait_cols, drop = FALSE]
    names(temp1) <- name_vec
    temp1 <- cbind(crosses_df, temp1)

    temp2 <- cv[, idx_cols, drop = FALSE]
    names(temp2) <- c(
      "GEBV.IDX",
      "TGV.IDX",
      "var.A.IDX",
      "SPV.IDX",
      "var.D.IDX",
      "var.TGV.IDX",
      "TSPV.IDX"
    )
    temp2 <- cbind(crosses_df, temp2)

    return(list(
      cross.df = temp1,
      index.df = temp2
    ))
  }

  out <- as.data.frame(temp)
  out <- out[, 1:(7 * ntraits), drop = FALSE]
  names(out) <- name_vec
  out <- cbind(crosses_df, out)
  return(out)
}
