#' Calculate optimal haploid values for candidate crosses
#'
#' Calculates the Optimal Haploid Value (OHV) of proposed biparental or
#' four-parent crosses from additive marker effects.
#'
#' The OHV describes the maximum additive genetic value obtainable by combining
#' favourable parental haplotypes and was proposed as a criterion for selecting
#' breeding material according to its potential to generate superior fixed
#' progeny (Daetwyler et al., 2015). Related haplotype-based approaches have
#' also been used to characterise Mendelian sampling effects and upper
#' selection limits (Cole and VanRaden, 2011).
#'
#' For each haplotype block, the function determines the most favourable
#' parental contribution according to the supplied marker effects and sums
#' these contributions across blocks. If \code{haplotype.blocks = NULL}, each
#' marker is treated as an independent block. User-defined haplotype blocks can
#' instead be supplied to account for sets of markers that should be inherited
#' jointly.
#'
#' The OHV represents an upper bound on additive cross potential rather than
#' the expected mean performance of progeny. In contrast to criteria such as
#' the superior progeny value (SPV), it does not use the expected progeny mean
#' and segregation variance to describe the distribution of progeny values.
#' Instead, it evaluates the most favourable combination of parental genetic
#' material that can be assembled under the specified haplotype-block
#' definition.
#'
#' Multiple traits can be evaluated simultaneously. If trait weights are
#' supplied, trait-specific OHVs are combined into a weighted OHV index.
#'
#' A practical consideration is that certain plant-breeding procedures based
#' on the OHV criterion have been the subject of patent protection, including
#' US11744199B2 (Daetwyler et al., 2023). Patent scope and status may differ
#' among jurisdictions and over time; users intending commercial application
#' should assess whether relevant patent claims apply to their intended use.
#'
#' @param crosses Matrix or data frame with two columns for biparental crosses
#'   or four columns for four-parent crosses. Parent identifiers may be integer
#'   row indices of \code{marker.mat} or character identifiers matching
#'   \code{rownames(marker.mat)}.
#' @param marker.mat Numeric marker matrix with individuals in rows and markers
#'   in columns.
#' @param marker.effects Numeric matrix of additive marker effects with markers
#'   in rows and traits in columns. Its number of rows must equal
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
#' @references
#' Cole, J. B. and VanRaden, P. M. (2011).
#' Use of haplotypes to estimate Mendelian sampling effects and selection
#' limits.
#' \emph{Journal of Animal Breeding and Genetics}, 128(6), 446--455.
#' \doi{10.1111/j.1439-0388.2011.00922.x}
#'
#' Daetwyler, H. D., Hayden, M. J., Spangenberg, G. C. and Hayes, B. J. (2015).
#' Selection on Optimal Haploid Value increases genetic gain and preserves more
#' genetic diversity relative to genomic selection.
#' \emph{Genetics}, 200(4), 1341--1348.
#' \doi{10.1534/genetics.115.178038}
#'
#' Daetwyler, H. D., Hayes, B. J., Robbins, K., Hayden, M. J. and
#' Spangenberg, G. (2023).
#' Selection based on optimal haploid value to create elite lines.
#' U.S. Patent US11744199B2.
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
