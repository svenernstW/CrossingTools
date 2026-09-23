#' Calculate segregation variance and superior progeny values for inbred lines
#'
#' Calculates expected genomic breeding values, segregation variances, and
#' superior progeny values for doubled haploid or recombinant inbred line
#' families derived from proposed two-way or four-way crosses.
#'
#' The calculation can account for additional generations of random mating
#' or inbreeding before line development. Multiple traits and optional weighted
#' trait indices are supported.
#'
#' Segregation (co)variances are calculated analytically from parental marker
#' genotypes, recombination fractions, and average allele-substitution effects.
#' For biparental crosses, the available methods follow the analytical
#' approaches of Lehermeier et al. (2017) and Osthushenrich et al. (2017).
#' Both provide expressions for doubled haploid (DH) and recombinant
#' inbred line (RIL) populations after an arbitrary number of intermating
#' generations. Four-way crosses are calculated using the extension to
#' multi-parental crosses described by Allier et al. (2019).
#'
#' Additive genetic values are expressed as breeding values using centred
#' marker genotypes, \eqn{M - 2p}, and average allele-substitution effects
#' (\eqn{\alpha}). The reference allele frequencies \eqn{p} can be supplied by
#' the user. If they are not supplied, they are calculated from
#' \code{marker.mat}. Thus, \code{p} defines the reference population relative
#' to which breeding values are expressed.
#'
#' Superior progeny values (SPV) combine the expected family breeding value
#' with the predicted segregation standard deviation and a user-supplied
#' selection intensity. For multiple traits, segregation covariance matrices
#' can additionally be calculated and combined with trait weights to obtain
#' an index-level SPV.
#'
#' @param crosses Matrix or data frame with two columns for two-way crosses or
#'   four columns for four-way crosses. Parent identifiers may be row indices of
#'   \code{marker.mat} or character identifiers matching
#'   \code{rownames(marker.mat)}.
#' @param genetic.map Data frame describing the genetic map. It must contain:
#'   \describe{
#'     \item{\code{site}}{Marker index or marker name matching
#'     \code{colnames(marker.mat)}.}
#'     \item{\code{chr}}{Chromosome identifier.}
#'     \item{\code{pos}}{Marker position on the chromosome in Morgan.}
#'   }
#'   All markers in \code{marker.mat} must be represented.
#' @param marker.mat Numeric marker dosage matrix with individuals in rows and
#'   markers in columns, coded 0, 1, and 2 for the counted allele.
#' @param marker.effects Numeric matrix of average allele-substitution effects
#'   (\eqn{\alpha}) with markers in rows and traits in columns. Effects should be
#'   parameterised relative to the reference population defined by \code{p}.
#'   If \code{p = NULL}, the reference allele frequencies are derived from
#'   \code{marker.mat}. Its number of rows must equal
#'   \code{ncol(marker.mat)}.
#' @param t Integer. Number of additional generations considered before final
#'   line development. Its interpretation depends on the selected offspring
#'   type and segregation-variance method.
#' @param intensity Numeric scalar giving the standardised selection intensity
#'   used to calculate superior progeny values. The default is 1.
#' @param type Character. Offspring type, either \code{"DH"} for doubled
#'   haploid lines or \code{"RIL"} for recombinant inbred lines.
#' @param weights Optional numeric vector with one weight per trait. When
#'   supplied, weighted index values are calculated for each cross.
#' @param p Optional numeric vector of reference allele frequencies, with one
#'   value per marker. These frequencies define the reference population used
#'   to centre marker genotypes for breeding-value calculations. If
#'   \code{NULL}, allele frequencies are calculated from \code{marker.mat}.
#' @param covariance Logical. If \code{TRUE}, also calculate segregation
#'   covariance matrices among traits for each cross.
#' @param method Method used to calculate segregation variance. Either
#'   \code{1} or \code{"lehermeier"} for the Lehermeier et al. (2017) method,
#'   or \code{2} or \code{"osthushenrich"} for the Osthushenrich et al. (2017)
#'   method. The latter is currently available only for two-way crosses.
#'   Four-way crosses are calculated using the multi-parental formulation of
#'   Allier et al. (2019).
#' @param nthreads Positive integer. Number of computational threads.
#'
#' @return If neither \code{weights} nor trait covariances are requested, a data
#'   frame containing the parental identifiers and, for each trait,
#'   \code{GEBV.<trait>}, \code{var.<trait>}, and \code{SPV.<trait>}.
#'
#'   If \code{weights} is supplied, a list additionally containing
#'   \code{index.df}, with columns \code{GEBV.IDX}, \code{var.IDX}, and
#'   \code{SPV.IDX}.
#'
#'   If \code{covariance = TRUE}, a list containing:
#'   \describe{
#'     \item{\code{cross.df}}{Cross-specific trait values and segregation
#'     variances.}
#'     \item{\code{index.df}}{Weighted index values, if \code{weights} is
#'     supplied.}
#'     \item{\code{covariances}}{Segregation covariance matrices for each
#'     cross.}
#'   }
#'
#' @references
#' Lehermeier, C., Teyssedre, S. and Schon, C.-C. (2017).
#' Genetic gain increases by applying the usefulness criterion with improved
#' variance prediction in selection of crosses.
#' \emph{Genetics}, 207(4), 1651--1661.
#' \doi{10.1534/genetics.117.300403}
#'
#' Osthushenrich, T., Frisch, M. and Herzog, E. (2017).
#' Genomic selection of crossing partners on basis of the expected mean and
#' variance of their derived lines.
#' \emph{PLOS ONE}, 12(12), e0188839.
#' \doi{10.1371/journal.pone.0188839}
#'
#' Allier, A., Moreau, L., Charcosset, A., Teyssedre, S. and
#' Lehermeier, C. (2019).
#' Usefulness criterion and post-selection parental contributions in
#' multi-parental crosses: Application to polygenic trait introgression.
#' \emph{G3: Genes|Genomes|Genetics}, 9(5), 1469--1479.
#' \doi{10.1534/g3.119.400129}
#'
#' @export

calc_spv_inbred <- function(crosses, genetic.map, marker.mat, marker.effects, t, intensity=1, type = "DH",
                            weights = NULL,p=NULL, covariance = FALSE,
                                   method = 1, nthreads = 4L) {

  marker.effects <- as.matrix(marker.effects)
  effects <- marker.effects

  traits <- colnames(marker.effects)

  if (is.null(traits)) {
    traits <- paste0("trait", seq_len(ncol(marker.effects)))
  }

  effects <- as.matrix(marker.effects)

  if (!is.numeric(intensity) ||
      length(intensity) != 1L ||
      !is.finite(intensity)) {
    stop("`intensity` must be a single finite numeric value.")
  }

  n.Threads <- nthreads
  if(!ncol(crosses) %in% c(2,4)){stop("ncol(crosses) needs to be 2 for two way crosses or 4 for three or four way crosses")}
  crosses_in <- crosses
  if (!isTRUE(all.equal(t, as.integer(t)))) {
    stop("t must be an integer.")
  }
  crosses <- as.matrix(crosses)

  if (is.numeric(crosses) || is.integer(crosses)) {
    if (any(!is.finite(crosses))) stop("`crosses` contains non-finite entries.")
    if (any(crosses < 1 | crosses > nrow(marker.mat))) {
      stop("Some genotype indices in `crosses` are outside 1..nrow(marker.mat).")
    }
    crosses2 <- as.data.frame(crosses)
  } else {
    if (is.null(rownames(marker.mat))) stop("Character `crosses` requires rownames(marker.mat).")

    idx <- match(as.vector(crosses), rownames(marker.mat))
    if (anyNA(idx)) stop("Some entries in `crosses` are not in rownames(marker.mat).")

    idx <- matrix(idx, nrow = nrow(crosses), ncol = ncol(crosses), byrow = FALSE)
    crosses2 <- as.data.frame(idx)
  }



  if (ncol(crosses) == 2){
    cross.type <- "2W"
  }  else {
    cross.type <- "4W"}


  if (cross.type == "2W") {
    names(crosses2) <- c("Parent1","Parent2")
  } else {
    names(crosses2) <- c("Parent1","Parent2","Parent3","Parent4")
    }

  crosses2 <- as.matrix(crosses2)
  storage.mode(crosses2) <- "integer"
  if (anyNA(crosses2)) stop("Internal error: `crosses2` contains NA after conversion.")

  if (is.numeric(method)) {
    if (method == 1) method <- "lehermeier"
    else if (method == 2) method <- "osthushenrich"
    else stop("`method` must be 1 or 2.")
  } else {
    method <- match.arg(method, c("lehermeier","osthushenrich"))
  }

  if(!type %in% c("DH","RIL")){stop("type needs to be one of DH or RIL!")}


  # ---- Normalize ----
  if (!is.matrix(marker.mat)) marker.mat <- as.matrix(marker.mat)
  if (!is.matrix(effects)) effects <- as.matrix(effects)
  crosses2 <- as.matrix(crosses2)
  ###############################################################################
  # Reference allele frequencies
  ###############################################################################

  if (is.null(p)) {

    # Default: the supplied marker population defines the reference population
    p <- colMeans(marker.mat) / 2

  } else {

    if (!is.numeric(p)) {
      stop("`p` must be a numeric vector of reference allele frequencies.")
    }

    p <- as.numeric(p)

    if (length(p) != ncol(marker.mat)) {
      stop(
        "`p` must contain one allele frequency per marker: ",
        "length(p) = ", length(p),
        ", ncol(marker.mat) = ", ncol(marker.mat), "."
      )
    }
  }

  if (any(!is.finite(p))) {
    stop("`p` must contain only finite values.")
  }

  if (any(p < 0 | p > 1)) {
    stop("All values in `p` must be between 0 and 1.")
  }

  ntraits <- ncol(marker.effects)
  if (ntraits == 1L) covariance <- FALSE

  if (ntraits == 1L && !is.null(weights)) {
    warning("Single trait: weights are ignored (index equals trait).")
    weights <- NULL
  }


  calculate.index <- !is.null(weights) && covariance
  calculate.simple.index <- !is.null(weights) && !covariance

  if (is.null(weights)) {
    weights <- rep(1, ncol(effects))
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != ncol(effects)) {
      stop("`weights` must have length equal to ncol(marker.effects).", call. = FALSE)
    }
    if (any(!is.finite(weights))) {
      stop("`weights` must contain only finite values.", call. = FALSE)
    }
  }

  if (calculate.simple.index) {
    effects.index <- effects %*% weights
    effects <- cbind(effects, effects.index)
  }




  req_cols <- c("site", "chr", "pos")
  if (!all(req_cols %in% names(genetic.map))){
    stop("`genetic.map` must contain columns: site, chr, pos.")
  }
  site <- genetic.map$site

  # Only needed if site is character
  if (is.numeric(site) || is.integer(site)) {

    site_idx <- as.integer(site)

    if (any(site_idx < 1 | site_idx > ncol(marker.mat), na.rm = TRUE))
      stop("`genetic.map$site` contains indices outside 1..ncol(marker.mat).")

    if (!all(seq_len(ncol(marker.mat)) %in% site_idx))
      stop("Some markers in marker.mat are missing from `genetic.map$site`.")

  } else {

    mnames <- colnames(marker.mat)
    if (is.null(mnames))
      stop("Character `genetic.map$site` requires colnames(marker.mat).")

    site_chr <- as.character(site)

    if (anyDuplicated(site_chr))
      stop("`genetic.map$site` contains duplicated marker names.")

    if (!all(mnames %in% site_chr))
      stop("Some markers in marker.mat are missing from `genetic.map$site`.")

    site_idx <- match(site_chr, mnames)
    if (anyNA(site_idx))
      stop("`genetic.map$site` contains names not found in colnames(marker.mat).")
  }

  # Order map by chr/pos and carry indices along
  o    <- order(genetic.map$chr, genetic.map$pos)
  map2 <- genetic.map[o, , drop = FALSE]
  ord  <- site_idx[o]

  marker.mat <- marker.mat[, ord, drop = FALSE]
  effects    <- effects[ord, , drop = FALSE]
  p          <- p[ord]

  # Relabel sites sequentially for downstream code
  map2$site <- seq_len(ncol(marker.mat))


  chr_levels <- unique(map2$chr)
  genmap_list <- lapply(chr_levels, function(cc) {
    as.matrix(map2$pos[map2$chr == cc])
  })


  if (length(n.Threads) != 1L || !is.finite(n.Threads) || n.Threads < 1 || n.Threads != as.integer(n.Threads)) {
    stop("`n.Threads` must be a positive integer.")
  }
  nThreads <- as.integer(n.Threads)

  if(type=="DH"){
    if(cross.type=="2W"){
      if (t < 0) {
        message("For two way crosses t needs to be >= 0, setting t = 0")
        t <- 0
      }
      if (method == "lehermeier") {
        temp <- cpp_calculate_covariance_lehermeier(
          Crosses    = crosses2,
          genMap     = genmap_list,
          M = marker.mat,
          U    = effects,
          t          = as.integer(t),
          intensity  = intensity,
          weights    = weights,
          p = p,
          covariance = covariance,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
      } else {
        temp <- cpp_calculate_covariance_osthushenrich(
          Crosses    = crosses2,
          genMap     = genmap_list,
          M = marker.mat,
          U    = effects,
          t          = as.integer(t),
          intensity  = intensity,
          weights    = weights,
          p = p,
          covariance = covariance,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
      }
    }

    if(cross.type=="4W"){
      if (t < 1) {
        message("For four way crosses t needs to be >= 1, setting t = 1")
        t <- 1
      }
      if (method == "lehermeier") {
        temp <- cpp_calculate_covariance_allier(
          Crosses    = crosses2,
          genMap     = genmap_list,
          M = marker.mat,
          U    = effects,
          t          = as.integer(t),
          intensity  = intensity,
          weights    = weights,
          p = p,
          covariance = covariance,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
      } else {
        print("Method 2 currently not supported for three and four way crosses, reverting to method = 1")
        temp <- cpp_calculate_covariance_allier(
          Crosses    = crosses2,
          genMap     = genmap_list,
          M = marker.mat,
          U    = effects,
          t          = as.integer(t),
          intensity  = intensity,
          weights    = weights,
          p = p,
          covariance = covariance,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
      }
    }
  }


  if(type=="RIL"){
    if (t < 0) {
      message("t needs to be >= 0, setting t = 0")
      t <- 0
    }
    if(cross.type=="2W"){
      if (method == "lehermeier") {
        temp <- cpp_calculate_covariance_RIL_lehermeier(
          Crosses    = crosses2,
          genMap     = genmap_list,
          M = marker.mat,
          U    = effects,
          t          = as.integer(t),
          intensity  = intensity,
          weights    = weights,
          p = p,
          covariance = covariance,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
      } else {
        temp <- cpp_calculate_covariance_RIL_osthushenrich(
          Crosses    = crosses2,
          genMap     = genmap_list,
          M = marker.mat,
          U    = effects,
          t          = as.integer(t),
          intensity  = intensity,
          weights    = weights,
          p = p,
          covariance = covariance,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
      }
    }

    if(cross.type=="4W"){
      if (t < 1) {
        message("For four way crosses t needs to be >= 1, setting t = 1")
        t <- 1
      }
      if (method == "lehermeier") {
        temp <- cpp_calculate_covariance_RIL_allier(
          Crosses    = crosses2,
          genMap     = genmap_list,
          M = marker.mat,
          U    = effects,
          t          = as.integer(t),
          intensity  = intensity,
          weights    = weights,
          p = p,
          covariance = covariance,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
      } else {
        print("Method 2 currently not supported for three and four way crosses, reverting to method = 1")
        temp <- cpp_calculate_covariance_RIL_allier(
          Crosses    = crosses2,
          genMap     = genmap_list,
          M = marker.mat,
          U    = effects,
          t          = as.integer(t),
          intensity  = intensity,
          weights    = weights,
          p = p,
          covariance = covariance,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
      }
    }
  }



  ntraits <- ncol(marker.effects)

  name_vec <- c(
    paste0("GEBV.",traits),
    paste0("var.",traits),
    paste0("SPV.",traits)
  )

  crosses_df <- as.data.frame(crosses_in, stringsAsFactors = FALSE)
  names(crosses_df) <- if (ncol(crosses) == 2) {
    c("parent1", "parent2")
  } else {
    c("parent1", "parent2", "parent3", "parent4")
  }

  if (covariance) {
    cv <- as.data.frame(temp$cross_values)

    # first 3*ntraits columns are always the trait-specific outputs
    temp1 <- cv[, 1:(3 * ntraits), drop = FALSE]
    names(temp1) <- name_vec
    temp1 <- cbind(crosses_df, temp1)

    if (calculate.index) {
      temp2 <- cv[, (3 * ntraits + 1):(3 * ntraits + 3), drop = FALSE]
      names(temp2) <- c("GEBV.IDX", "var.IDX", "SPV.IDX")
      temp2 <- cbind(crosses_df, temp2)

      if (isTRUE(temp$check_psd)) {
        warning("Some segregation covariance matrices were not psd; the nearest psd matrix was used for index calculation.")
      }

      return(list(
        cross.df = temp1,
        index.df = temp2,
        covariances = temp$covariances
      ))
    } else {
      return(list(
        cross.df = temp1,
        covariances = temp$covariances
      ))
    }
  }

  if (calculate.simple.index) {
    cv <- as.data.frame(temp)

    ntraits.cpp <- ncol(effects)  # original traits + 1 appended index trait

    trait_cols <- c(
      seq_len(ntraits),
      ntraits.cpp + seq_len(ntraits),
      2 * ntraits.cpp + seq_len(ntraits)
    )

    idx_cols <- c(
      ntraits.cpp,
      2 * ntraits.cpp,
      3 * ntraits.cpp
    )

    temp1 <- cv[, trait_cols, drop = FALSE]
    names(temp1) <- name_vec
    temp1 <- cbind(crosses_df, temp1)

    temp2 <- cv[, idx_cols, drop = FALSE]
    names(temp2) <- c("GEBV.IDX", "var.IDX", "SPV.IDX")
    temp2 <- cbind(crosses_df, temp2)

    return(list(
      cross.df = temp1,
      index.df = temp2
    ))
  }

  out <- as.data.frame(temp)
  out <- out[, 1:(3 * ntraits), drop = FALSE]
  names(out) <- name_vec
  out <- cbind(crosses_df, out)
  return(out)
}
