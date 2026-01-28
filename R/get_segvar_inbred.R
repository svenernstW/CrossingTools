#' Segregation variance for inbred lines from specified crosses
#'
#' Calculate expected genomic estimated breeding values (EGEBV), segregation
#' variance, and superior progeny value (SPV) for doubled haploid or recombinant inbred lines
#' #' derived after \code{t} rounds of random mating. works for two way, three way and four way crosses
#'
#' @param crosses A matrix or data.frame (n_crosses x 2 for two way crosses or n_crosses x 4 for four way crosses) specifying the parents for each
#'   proposed cross. Entries may be either:
#'   \itemize{
#'     \item integer indices referring to rows of \code{marker.mat}, or
#'     \item character identifiers matching \code{rownames(marker.mat)}.
#'   }
#' @param marker.mat Numeric marker matrix (markers in columns, genotypes in rows) with
#'   entries in \code{c(0, 2)} corresponding to \code{c("AA", "BB")}.
#' @param genetic.map A data.frame with columns:
#'   \itemize{
#'     \item \code{site}: integer marker index (1..ncol(marker.mat))
#'     \item \code{chr}: chromosome identifier (numeric)
#'     \item \code{pos}: position on the chromosome in morgan
#'   }
#'   All markers in \code{marker.mat} must appear in \code{genetic.map$site}.
#' @param effects Numeric matrix of marker effects (markers in rows, traits in columns).
#'   Must have \code{nrow(effects) == ncol(marker.mat)}.
#' @param t Integer. Number of random‐mating generations before DH creation.
#' @param intensity Double. Standardized selection differential, efaults to 1.
#' @param covariance Logical. If \code{TRUE}, also compute the segregation
#'   covariance matrix between all supplied traits.
#' @param weights Numeric vector of length \code{ncol(effects)} with fixed trait weights.
#' If supplied the function calculates index values for each cross. Only works if \code{covariance = TRUE}; otherwise ignored.
#' @param method Character. Which method to use; one of \code{c(1, 2)} for
#'   Lehermeier or Osthushenrich.
#' @param type intended offspring typew can be "RIL" for recombinant inbred lines or "DH" for double haploid lines.
#' @param nthreads Integer (default 4). Number of OpenMP threads (if enabled at compile time).
#'
#' @return If \code{covariance = TRUE}, a list with element \code{cross.df}
#'   (a data.frame) whose columns are named \code{EGBEV1, var1, SPV1, EGEBV2, ...} and a list with element \code{covariances}.
#'   If \code{covariance = TRUE} and \code{calculate.index = TRUE}, additionally an element \code{index.df}
#'   (a data.frame) with \code{IDX, VAR.IDX, SPV.IDX}.
#'   If \code{covariance = FALSE}, a data.frame with the same columns \code{EGEBV1, var1, SPV1, EGBEV2, ...}.
#' @export
get_segvar_inbred <- function(crosses, genetic.map, marker.mat, effects, t, intensity=NULL, type = "DH",
                                   covariance = FALSE,  weights = NULL,
                                   method = 1, nthreads = 4L) {
  n.Threads <- nthreads
  if(!ncol(crosses) %in% c(2,4)){stop("ncol(crosses) needs to be 2 for two way crosses or 4 for three or four way crosses")}
  crosses_in <- crosses

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
  calculate.index <- !is.null(weights)

  if (ncol(effects) == 1 & is.null(weights)) { weights <- c(1) }

    # ---- Checks ----
  if (ncol(marker.mat) <= 0L) stop("marker.mat must have markers in columns.")
  if (nrow(effects) != ncol(marker.mat)) stop("effects must have nrow(effects) == ncol(marker.mat).")
  if (!is.logical(covariance) || length(covariance) != 1L) stop("`covariance` must be TRUE/FALSE.")
  if (!is.numeric(t) || length(t) != 1L || t < 0 || abs(t - round(t)) > .Machine$double.eps^0.5)
    stop("`t` must be a single non-negative integer.")
  if (!is.numeric(intensity) || length(intensity) != 1L){
    warning("`intensity` not provided, using 1.")
  intensity <- 1
  }
  if (calculate.index && is.null(weights)) {
    stop("`weights` is required when calculate.index = TRUE.")
  }
  if (calculate.index && length(weights) != ncol(effects)) {
    stop("`weights` must have length equal to ncol(effects) when calculate.index = TRUE.")
  }
  if (calculate.index) {
    weights <- as.numeric(weights)
    if (any(!is.finite(weights))) stop("`weights` must contain only finite values.")
  }

  if (t <= 0) {
    message("t needs to be >= 1, setting t = 1")
    t <- 1
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

  # relabel sites to 1..p for downstream code
  map2$site <- seq_len(ncol(marker.mat))

  chr_levels <- unique(map2$chr)
  genmap_list <- lapply(chr_levels, function(cc) {
    as.matrix(map2$pos[map2$chr == cc])
  })

  if (ncol(effects) == 1L) covariance <- FALSE

  if (!covariance && calculate.index) {
    message("covariance == FALSE: ignoring index related arguments")
    calculate.index <- FALSE
    weights <- rep(0, ncol(effects))
  }

  if (length(n.Threads) != 1L || !is.finite(n.Threads) || n.Threads < 1 || n.Threads != as.integer(n.Threads)) {
    stop("`n.Threads` must be a positive integer.")
  }
  nThreads <- as.integer(n.Threads)

  if(type=="DH"){
    if(cross.type=="2W"){
      if (method == "lehermeier") {
        temp <- cpp_calculate_covariance_lehermeier(
          Crosses    = crosses2,
          genMap     = genmap_list,
          M = marker.mat,
          U    = effects,
          t          = as.integer(t),
          intensity  = intensity,
          weights    = weights,
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
          covariance = covariance,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
      }
    }

    if(cross.type=="4W"){
      if (method == "lehermeier") {
        temp <- cpp_calculate_covariance_allier(
          Crosses    = crosses2,
          genMap     = genmap_list,
          M = marker.mat,
          U    = effects,
          t          = as.integer(t),
          intensity  = intensity,
          weights    = weights,
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
          covariance = covariance,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
      }
    }
  }


  if(type=="RIL"){
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
          covariance = covariance,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
      }
    }

    if(cross.type=="4W"){
      if (method == "lehermeier") {
        temp <- cpp_calculate_covariance_RIL_allier(
          Crosses    = crosses2,
          genMap     = genmap_list,
          M = marker.mat,
          U    = effects,
          t          = as.integer(t),
          intensity  = intensity,
          weights    = weights,
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
          covariance = covariance,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
      }
    }
  }



  name_vec <- paste0(rep(c("EGEBV","var","SPV"), each = ncol(effects)), seq_len(ncol(effects)))
  crosses_df <- as.data.frame(crosses_in, stringsAsFactors = FALSE)
  names(crosses_df) <- if (ncol(crosses) == 2) c("parent1","parent2") else c("parent1","parent2","parent3","parent4")


  if (covariance) {
    if (calculate.index) {
      temp1 <- as.data.frame(temp$cross_values)[, 1:(3 * ncol(effects)), drop = FALSE]
      names(temp1) <- name_vec
      temp1 <- cbind(crosses_df, temp1)

      temp2 <- as.data.frame(temp$cross_values)[, ((3 * ncol(effects)) + 1):ncol(as.data.frame(temp$cross_values)), drop = FALSE]
      names(temp2) <- c("IDX", "var.IDX", "SPV.IDX")
      temp2 <- cbind(crosses_df, temp2)
      if (temp$check_psd) {
        warning("Some segregation covariance matrices were not psd; for index calculation the nearest psd was used in these places")
      }
      return(list(cross.df = temp1, index.df = temp2, covariances = temp$covariances))
    } else {
      temp1 <- as.data.frame(temp$cross_values)[, 1:(3 * ncol(effects)), drop = FALSE]
      names(temp1) <- name_vec
      temp1 <- cbind(crosses_df, temp1)
      return(list(cross.df = temp1, covariances = temp$covariances))
    }
  } else {
    out <- as.data.frame(temp)
    names(out) <- name_vec
    out <- cbind(crosses_df, out)
    return(out)
  }
}
