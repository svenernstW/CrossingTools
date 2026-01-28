#' Optimal Haploid Value (OHV) for proposed crosses
#'
#' Computes the Optimal Haploid Value (OHV; Daetwyler et al.) for each proposed cross
#' using block-wise local genomic estimated breeding values (GEBVs).
#'
#' For each haplotype block, a local GEBV is calculated for every parent as the sum of
#' marker genotypes within the block weighted by their corresponding marker effects.
#' The OHV of a cross is then obtained by summing, across blocks, the maximum of the two
#' parents’ local block GEBVs.
#'
#' If \code{haplotype.blocks} is not supplied , each marker column of
#' \code{marker.mat} is treated as an independent block.
#'
#' @param crosses A matrix or data.frame (n_crosses x 2 for two way crosses or n_crosses x 4 for four way crosses) specifying the parents for each
#'   proposed cross. Entries may be either:
#'   \itemize{
#'     \item integer indices referring to rows of \code{marker.mat}, or
#'     \item character identifiers matching \code{rownames(marker.mat)}.
#'   }
#'
#' @param haplotype.blocks Optional data.frame defining haplotype blocks with columns:
#'   \itemize{
#'     \item \code{block}: block identifier (integer, character, or factor). Markers sharing
#'       the same value belong to the same haplotype block.
#'     \item \code{site}: marker identifier within \code{marker.mat}, given either as an
#'       integer column index (1..ncol(marker.mat)) or as a marker name matching
#'       \code{names(marker.mat)}.
#'   }
#'   Each marker may appear in at most one block.
#'
#' @param marker.mat Numeric genotype or marker matrix with individuals in rows and
#'   markers in columns.
#'
#' @param effects Numeric matrix of marker effects with
#'   \code{nrow(effects) == ncol(marker.mat)} with traits in columns.
#'
#' @param weights Optional numeric vector of length \code{ncol(effects)}. If supplied,
#'   a weighted index is computed as a linear combination of the trait-specific OHVs.
#'
#' @param nthreads Integer >= 1. Number of threads used for parallel computation.
#'
#' @return If \code{weights} is \code{NULL}, a data.frame containing the cross definition
#'   columns followed by \code{OHV.1, OHV.2, ...} for each trait.
#'
#'   If \code{weights} is supplied, a list with elements:
#'   \itemize{
#'     \item \code{cross.df}: data.frame of crosses and trait-specific OHVs
#'     \item \code{index.df}: data.frame of crosses and the weighted OHV index
#'   }
#'
#' @export

get_optimal_haploid_value <- function(crosses, haplotype.blocks = NULL,
                                            marker.mat, effects, weights = NULL, nthreads = 4L) {
  n.Threads  <- nthreads
  # ---- Normalize inputs ----
  if (!is.matrix(marker.mat)) marker.mat <- as.matrix(marker.mat)
  effects <- as.matrix(effects)
  crosses <- as.matrix(crosses)

  n.Threads <- nthreads
  if(!ncol(crosses) %in% c(2,4)){stop("ncol(crosses) needs to be 2 for two way crosses or 4 for three or four way crosses")}
  if (ncol(crosses) == 2){
    cross.type <- "2W"
  }  else {
    cross.type <- "4W"}


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






  if (cross.type == "2W") {
    names(crosses2) <- c("Parent1","Parent2")
  } else {
    names(crosses2) <- c("Parent1","Parent2","Parent3","Parent4")
  }

  crosses2 <- as.matrix(crosses2)
  storage.mode(crosses2) <- "integer"
  if (anyNA(crosses2)) stop("Internal error: `crosses2` contains NA after conversion.")


  # ---- Basic checks: marker.mat / effects ----
  if (!is.numeric(marker.mat)) stop("`marker.mat` must be a numeric matrix.")
  if (ncol(marker.mat) < 1L || nrow(marker.mat) < 1L) stop("`marker.mat` must have at least 1 row and 1 column.")

  if (!is.numeric(effects)) stop("`effects` must be numeric.")
  if (nrow(effects) != ncol(marker.mat)) {
    stop("`effects` must have nrow(effects) == ncol(marker.mat). Got nrow(effects) = ",
         nrow(effects), " but ncol(marker.mat) = ", ncol(marker.mat), ".")
  }
  if (any(!is.finite(effects))) stop("`effects` contains non-finite values.")
  if (ncol(effects) < 1L) stop("`effects` must have at least one trait column.")


  # ---- n.Threads ----
  if (length(n.Threads) != 1L || !is.finite(n.Threads) ||
      abs(n.Threads - round(n.Threads)) > .Machine$double.eps^0.5 || n.Threads < 1) {
    stop("`n.Threads` must be a single integer >= 1.")
  }
  nThreads <- as.integer(n.Threads)

  # ---- haplotype.blocks: data.frame(block, site) ----
  p <- ncol(marker.mat)
  if (is.null(haplotype.blocks) || (is.data.frame(haplotype.blocks) && nrow(haplotype.blocks) == 0L)) {
    HBlocks <- lapply(seq_len(p), function(j) as.integer(j))
  } else {
    if (!is.data.frame(haplotype.blocks)) {
      stop("`haplotype.blocks` must be a data.frame with columns `block` and `site` (or NULL).")
    }
    if (!all(c("block", "site") %in% names(haplotype.blocks))) {
      stop("`haplotype.blocks` must contain columns: `block` and `site`.")
    }
    if (nrow(haplotype.blocks) < 1L) stop("`haplotype.blocks` has 0 rows; use NULL or provide at least one row.")

    hb <- haplotype.blocks[, c("block", "site")]
    if (any(is.na(hb$block))) stop("`haplotype.blocks$block` contains NA.")
    if (any(is.na(hb$site)))  stop("`haplotype.blocks$site` contains NA.")

    # Convert site -> integer indices into marker.mat columns
    if (is.numeric(hb$site) || is.integer(hb$site)) {
      site_idx <- as.integer(hb$site)
      if (any(!is.finite(site_idx))) stop("`haplotype.blocks$site` contains non-finite values.")
      if (any(site_idx < 1 | site_idx > p)) stop("`haplotype.blocks$site` contains indices outside 1..ncol(marker.mat).")
    } else {
      mnames <- colnames(marker.mat)
      if (is.null(mnames)) stop("Character `haplotype.blocks$site` requires colnames(marker.mat).")
      site_chr <- as.character(hb$site)
      site_idx <- match(site_chr, mnames)
      if (anyNA(site_idx)) {
        bad <- unique(site_chr[is.na(site_idx)])
        stop("Some `haplotype.blocks$site` names are not in colnames(marker.mat): ",
             paste(bad, collapse = ", "))
      }
      site_idx <- as.integer(site_idx)
    }

    hb$site_idx <- site_idx

    # Each marker appears at most once
    if (any(duplicated(hb$site_idx))) {
      dup_sites <- unique(hb$site[duplicated(hb$site_idx)])
      stop("Some markers appear in multiple blocks: ", paste(dup_sites, collapse = ", "))
    }

    # Split into list by block, preserving first-appearance order
    ord_blocks <- unique(hb$block)
    HBlocks <- lapply(ord_blocks, function(b) as.integer(hb$site_idx[hb$block == b]))
  }

  if(cross.type=="2W"){
    # ---- Call C++ for each trait ----
    temp <- vector("list", ncol(effects))
    for (i in seq_len(ncol(effects))) {
      temp[[i]] <- cpp_calcOHV(
        Crosses  = crosses2,
        HBlocks  = HBlocks,
        M        = marker.mat,
        mu_vec   = effects[, i],
        nThreads = n.Threads
      )
    }

    out <- as.data.frame(do.call(cbind, temp))
    names(out) <- paste0("OHV.", seq_len(ncol(effects)))

  }

  if(cross.type=="4W"){
    # ---- Call C++ for each trait ----
    temp <- vector("list", ncol(effects))
    for (i in seq_len(ncol(effects))) {
      temp[[i]] <- cpp_calcOHV4W(
        Crosses  = crosses2,
        HBlocks  = HBlocks,
        M        = marker.mat,
        mu_vec   = effects[, i],
        nThreads = n.Threads
      )
    }

    out <- as.data.frame(do.call(cbind, temp))
    names(out) <- paste0("OHV.", seq_len(ncol(effects)))

  }

  if(!is.null(weights)){
    if(ncol(effects)!=length(weights)){
      stop("if weights are provided they have to have the same length as ncol(effects)")
    }
    idx <- out[[1]] * weights[1]
    for (j in 2:length(weights)) {
      idx <- idx + out[[j]] * weights[j]
    }
    idx <- data.frame(index = idx)

    out <- list(cross.df=cbind(as.data.frame(crosses_in, stringsAsFactors = FALSE), out),index.df=cbind(crosses,idx))
    out
  }else{
    out <- cbind(as.data.frame(crosses_in, stringsAsFactors = FALSE), out)
  }
  # keep original crosses provided by user (character or numeric) in output

  out
}
