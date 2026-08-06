#' Calculates the Optimal Haploid Value (OHV) of proposed two-way or four-way
#' crosses from marker effects. If haplotype blocks are provided, marker effects
#' are first aggregated within blocks; otherwise, each marker is treated as an
#' independent block. Trait-specific OHVs and an optional weighted OHV index are
#' returned.
#'
#' @param crosses Matrix or data frame with two columns for two-way crosses or
#'   four columns for four-way crosses. Parent identifiers may be
#'   integer row indices of \code{marker.mat} or character identifiers matching
#'   \code{rownames(marker.mat)}.
#' @param marker.mat Numeric marker matrix with individuals in rows and markers
#'   in columns.
#' @param marker.effects Numeric matrix of marker effects with markers in rows
#'   and traits in columns. Its number of rows must equal
#'   \code{ncol(marker.mat)}.
#' @param weights Optional numeric vector with one weight per trait. When
#'   supplied, a weighted OHV index is calculated as the linear combination of
#'   the trait-specific OHVs.
#' @param haplotype.blocks Optional data frame defining the haplotype blocks. It
#'   must contain:
#'   \describe{
#'     \item{\code{block}}{Block identifier. Markers with the same identifier
#'     belong to the same haplotype block.}
#'     \item{\code{site}}{Marker identifier given as a column index of
#'     \code{marker.mat} or as a marker name matching
#'     \code{colnames(marker.mat)}.}
#'   }
#'   Each marker may occur in at most one block. If \code{NULL}, each marker is
#'   treated as an independent block.
#' @param nthreads Positive integer. Number of computational threads.
#'
#' @return If \code{weights = NULL}, a data frame containing the original cross
#'   definitions followed by one \code{OHV.<trait>} column per trait.
#'
#'   If \code{weights} is supplied, a list containing:
#'   \describe{
#'     \item{\code{cross.df}}{The original cross definitions and the
#'     trait-specific OHVs.}
#'     \item{\code{index.df}}{The cross definitions and the weighted OHV index
#'     in column \code{OHV.IDX}.}
#'   }
#'
#' @export

calc_optimal_haploid_value <- function(crosses,
                                            marker.mat, marker.effects, weights = NULL,haplotype.blocks = NULL, nthreads = 4L) {

  traits <- names(marker.effects)
  n.Threads  <- nthreads
  effects <- marker.effects
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
    names(out) <- paste0("OHV.", traits)

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
    names(out) <- paste0("OHV.", traits)

  }

  if(!is.null(weights)){
    if(ncol(effects)!=length(weights)){
      stop("if weights are provided they have to have the same length as ncol(effects)")
    }
    idx <- out[[1]] * weights[1]
    for (j in 2:length(weights)) {
      idx <- idx + out[[j]] * weights[j]
    }
    idx <- data.frame(OHV.IDX = idx)

    out <- list(cross.df=cbind(as.data.frame(crosses_in, stringsAsFactors = FALSE), out),index.df=cbind(crosses,idx))
    out
  }else{
    out <- cbind(as.data.frame(crosses_in, stringsAsFactors = FALSE), out)
  }
  # keep original crosses provided by user (character or numeric) in output

  out
}
