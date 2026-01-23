#' Segregation variance for DH lines from specified four- or three-way crosses (Allier)
#'
#' Calculate expected genomic estimated breeding values (EGBV), segregation
#' variance, and superior progeny value (SPV) for doubled haploid (DH) lines
#' derived after \code{t} rounds of random mating.
#'
#' @param crosses A data.frame or matrix (n_crosses x 4) of parental indices
#'   (row indices into \code{marker.mat}) specifying four- or three-way crosses.
#'   For three-way crosses, the first two crossing partners are the same.
#'   You can also compute the variance for two heterozygous individuals by
#'   passing haplotypes in \code{marker.mat}; then indices refer to haplotypes.
#' @param marker.mat Numeric marker matrix (genotypes/haplotypes in rows, markers in columns)
#'   with entries in \code{c(0, 2)} corresponding to \code{c("AA", "BB")}.
#'   For heterozygotes, store haplotypes in \code{marker.mat} (each individual uses two rows).
#' @param genetic.map A data.frame with columns:
#'   \itemize{
#'     \item \code{site}: integer marker index (1..ncol(marker.mat))
#'     \item \code{chr}: chromosome identifier
#'     \item \code{pos}: position on the chromosome in morgan
#'   }
#'   All markers in \code{marker.mat} must appear in \code{genetic.map$site}.
#' @param effects Numeric matrix of marker effects (markers in rows, traits in columns).
#'   Must have \code{nrow(effects) == ncol(marker.mat)}.
#' @param t Integer. Number of random-mating generations before DH creation.
#' @param intensity Numeric. Standardized selection differential.
#' @param covariance Logical. If \code{TRUE}, also compute the segregation
#'   covariance matrix between all supplied traits.
#' @param weights Numeric vector of length \code{ncol(effects)} with fixed trait weights.
#' If supplied the function calculates index values for each cross. Only works if \code{covariance = TRUE}; otherwise ignored.
#' @param n.Threads Integer (default 4). Number of OpenMP threads (if enabled at compile time).
#'
#' @return If \code{covariance = TRUE}, a list with element \code{cross_values}
#'   (a data.frame) whose columns are named \code{EGBV1, var1, SPV1, EGBV2, ...} and a list with element \code{covariances} with segregation covariance matrix for each cross.
#'   If \code{covariance = TRUE}, a list with element \code{cross_values}
#'   (a data.frame) whose columns are named \code{EGBV1, var1, SPV1, EGBV2, ...}, element \code{index} (a data.frame) with \code{EGBV_DG, varDG, SPVDG} for the  index
#'   and a list with element \code{covariances} with segregation covariance matrix for each cross.
#'   If \code{covariance = FALSE}, a data.frame with the same columns \code{EGBV1, var1, SPV1, EGBV2, ...}.
#'
#' @examples
#' \dontrun{
#' out <- calculate_variances_DH(
#'   crosses = matrix(c(1,2, 3,4), ncol = 4, byrow = TRUE),
#'   genetic.map = data.frame(site = 1:ncol(marker.mat), chr = 1, pos = seq_len(ncol(marker.mat))),
#'   marker.mat = marker.mat,
#'   effects = matrix(rnorm(ncol(marker.mat)), ncol = 1),
#'   t = 0,
#'   intensity = 1.0,
#'   covariance = FALSE
#' )
#' }
#' @export
calculate_variances_4W_DH <- function(crosses, genetic.map, marker.mat, effects, t, intensity,
                                      covariance = FALSE, weights = NULL,
                                      n.Threads = 4L) {
  # ---- Normalize ----
  if (!is.matrix(marker.mat)) marker.mat <- as.matrix(marker.mat)
  if (!is.matrix(effects)) effects <- as.matrix(effects)
  crosses <- as.matrix(crosses)
  calculate.index <- !is.null(weights)

  if (ncol(effects) == 1 & is.null(weights)) { weights <- c(1) }
  # ---- Checks ----
  if (ncol(marker.mat) <= 0L) stop("marker.mat must have markers in columns.")
  if (nrow(effects) != ncol(marker.mat)) stop("effects must have nrow(effects) == ncol(marker.mat).")
  if (!is.logical(covariance) || length(covariance) != 1L) stop("`covariance` must be TRUE/FALSE.")
  if (!is.numeric(t) || length(t) != 1L || t < 0 || abs(t - round(t)) > .Machine$double.eps^0.5)
    stop("`t` must be a single non-negative integer.")
  if (!is.numeric(intensity) || length(intensity) != 1L)
    stop("`intensity` must be a single numeric.")
  if (ncol(crosses) != 4L) stop("`crosses` must have 4 columns (four-/three-way).")

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
    print("t needs to be larger or equal to 1, seting t=1")
    t <- 1
  }

  genos <- unique(as.vector(crosses))
  if (any(!is.finite(genos))) stop("`crosses` contains non-finite entries.")
  if (any(genos < 1 | genos > nrow(marker.mat))) stop("Some genotype indices in `crosses` are outside 1..nrow(marker.mat).")

  req_cols <- c("site", "chr", "pos")
  if (!all(req_cols %in% names(genetic.map)))
    stop("`genetic.map` must contain columns: site, chr, pos.")
  if (any(genetic.map$site < 1 | genetic.map$site > ncol(marker.mat)))
    stop("`genetic.map$site` contains indices outside 1..ncol(marker.mat).")
  if (!all(seq_len(ncol(marker.mat)) %in% genetic.map$site))
    stop("Some markers in marker.mat are missing from `genetic.map$site`.")

  map2 <- genetic.map[order(genetic.map$chr, genetic.map$pos), , drop = FALSE]
  ord  <- as.integer(map2$site)

  marker.mat <- marker.mat[, ord, drop = FALSE]
  effects <- effects[ord,  , drop = FALSE]

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

  temp <- cpp_calculate_covariance_allier(
    Crosses   = crosses,
    genMap    = genmap_list,
    M         = marker.mat,
    U         = effects,
    t         = as.numeric(t),
    intensity = intensity,
    covariance = covariance,
    calcindex = calculate.index,
    weights = weights,
    nThreads  = nThreads
  )

  name_vec <- paste0(rep(c("EGBV","var","SPV"), each = ncol(effects)), seq_len(ncol(effects)))
  crosses <- as.data.frame(crosses)
  names(crosses) <- c("parent1","parent2","parent3","parent4")

  if (covariance) {
    if (calculate.index) {
      temp1 <- as.data.frame(temp$cross_values)[, 1:(3 * ncol(effects))]
      names(temp1) <- name_vec
      temp1 <- cbind(crosses, temp1)
      temp2 <- as.data.frame(temp$cross_values)[, ((3 * ncol(effects)) + 1):ncol(as.data.frame(temp$cross_values))]
      names(temp2) <- c("IDG_A","VARIDG_A","SPVIDG_A")
      temp2 <- cbind(crosses, temp2)
      if (temp$check_psd) {
        warning("Some segregation covariance matrices were not psd; for index calculation the nearest psd was used in these places")
      }
      return(list(cross_values = temp1, index = temp2, covariances = temp$covariances))
    } else {
      temp1 <- as.data.frame(temp$cross_values)[, 1:(3 * ncol(effects))]
      names(temp1) <- name_vec
      temp1 <- cbind(crosses, temp1)
      return(list(cross_values = temp1, covariances = temp$covariances))
    }
  } else {
    out <- as.data.frame(temp)
    names(out) <- name_vec
    out <- cbind(crosses, out)
    return(out)
  }
}
