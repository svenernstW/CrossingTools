#' Segregation variance for RIL  from specified four- or three-way crosses (Allier)
#'
#' Calculate expected genomic estimated breeding values (EGBV), segregation
#' variance, and superior progeny value (SPV) for  recombinant inbred lines(RIL)
#' derived after \code{t} rounds of random mating.
#'
#' @param crosses A data.frame or matrix (n_crosses x 4) of parental indices
#'   (row indices into \code{M}) specifying four- or three-way crosses.
#'   For three-way crosses, the first two crossing partners are the same.
#'   You can also compute the variance for two heterozygous individuals by
#'   passing haplotypes in \code{M}; then indices refer to haplotypes.
#' @param M Numeric marker matrix (genotypes/haplotypes in rows, markers in columns)
#'   with entries in \code{c(0, 2)} corresponding to \code{c("AA", "BB")}.
#'   For heterozygotes, store haplotypes in \code{M} (each individual uses two rows).
#' @param genetic.map A data.frame with columns:
#'   \itemize{
#'     \item \code{site}: integer marker index (1..ncol(M))
#'     \item \code{chr}: chromosome identifier
#'     \item \code{pos}: position on the chromosome in cM
#'   }
#'   All markers in \code{M} must appear in \code{genetic.map$site}.
#' @param U Numeric matrix of marker effects (markers in rows, traits in columns).
#'   Must have \code{nrow(U) == ncol(M)}.
#' @param t Integer. Number of random-mating generations before RIL creation.
#' @param intensity Numeric. Standardized selection differential.
#' @param covariance Logical. If \code{TRUE}, also compute the segregation
#'   covariance matrix between all supplied traits.
#' @param calculate.gains Logical. If \code{TRUE}, calculate the desired gains index
#'  based on segregation (co)variances. This only works if covariance is TRUE, else will be ignored.
#' @param gains a vector of length equal to the number of traits with values representing the desired gains
#'
#' @param n.Threads Integer (default 4). Number of OpenMP threads (if enabled at compile time).
#'
#' @return If \code{covariance = TRUE}, a list with element \code{cross_values}
#'   (a data.frame) whose columns are named \code{EGBV1, var1, SPV1, EGBV2, ...} and a list with element \code{covariances} with segregation covariance matrix for each cross.
#'   If \code{covariance = TRUE}, a list with element \code{cross_values}
#'   (a data.frame) whose columns are named \code{EGBV1, var1, SPV1, EGBV2, ...}, element \code{gains} (a data.frame) with \code{EGBV_DG, varDG, SPVDG} for the desired gains index
#'   and a list with element \code{covariances} with segregation covariance matrix for each cross.
#'   If \code{covariance = FALSE}, a data.frame with the same columns \code{EGBV1, var1, SPV1, EGBV2, ...}.
#'
#' @examples
#' \dontrun{
#' out <- calculate_variances_RIL(
#'   crosses = matrix(c(1,2, 3,4), ncol = 4, byrow = TRUE),
#'   genetic.map = data.frame(site = 1:ncol(M), chr = 1, pos = seq_len(ncol(M))),
#'   M = M,
#'   U = matrix(rnorm(ncol(M)), ncol = 1),
#'   t = 0,
#'   intensity = 1.0,
#'   covariance = FALSE
#' )
#' }
#' @export
calculate_variances_4W_RIL <- function(crosses, genetic.map, M, U, t, intensity,
                                       covariance = FALSE, calculate.gains = FALSE, gains = NULL,
                                       n.Threads = 4L) {
  # ---- Normalize ----
  if (!is.matrix(M)) M <- as.matrix(M)
  if (!is.matrix(U)) U <- as.matrix(U)
  crosses <- as.matrix(crosses)

  # ---- Checks ----
  if (ncol(M) <= 0L) stop("M must have markers in columns.")
  if (nrow(U) != ncol(M)) stop("U must have nrow(U) == ncol(M).")
  if (!is.logical(covariance) || length(covariance) != 1L) stop("`covariance` must be TRUE/FALSE.")
  if (!is.numeric(t) || length(t) != 1L || t < 0 || abs(t - round(t)) > .Machine$double.eps^0.5)
    stop("`t` must be a single non-negative integer.")
  if (!is.numeric(intensity) || length(intensity) != 1L)
    stop("`intensity` must be a single numeric.")
  if (ncol(crosses) != 4L) stop("`crosses` must have 4 columns (four-/three-way).")
  if(!covariance & calculate.gains){
    print("covariance == FALSE, ignoring desired gains related arguments")
    calculate.gains <- FALSE
    gains <- rep(0,ncol(U))
  }

  if(calculate.gains & is.null((gains))){
    print("gains == NULL, ignoring desired gains related arguments")
    calculate.gains <- FALSE
    gains <- rep(0,ncol(U))
  }

  if(length(gains)!= ncol(U)){
    print("length(gains)!= ncol(U), ignoring desired gains related arguments")
    calculate.gains <- FALSE
    gains <- rep(0,ncol(U))
  }

  if(t<=0){
    print("t needs to be larger or equal to 1, seting t=1")
    t=1
  }

  # Validate crosses
  genos <- unique(as.vector(crosses))
  if (any(!is.finite(genos))) stop("`crosses` contains non-finite entries.")
  if (any(genos < 1 | genos > nrow(M))) stop("Some genotype indices in `crosses` are outside 1..nrow(M).")

  # Validate map
  req_cols <- c("site", "chr", "pos")
  if (!all(req_cols %in% names(genetic.map)))
    stop("`genetic.map` must contain columns: site, chr, pos.")
  if (any(genetic.map$site < 1 | genetic.map$site > ncol(M)))
    stop("`genetic.map$site` contains indices outside 1..ncol(M).")
  if (!all(seq_len(ncol(M)) %in% genetic.map$site))
    stop("Some markers in M are missing from `genetic.map$site`.")

  # ---- Crucial: reorder M *and* U to (chr, pos), and rebuild map consistently ----
  map2 <- genetic.map[order(genetic.map$chr, genetic.map$pos), , drop = FALSE]
  ord  <- as.integer(map2$site)  # these are original column indices

  M <- M[, ord, drop = FALSE]
  U <- U[ord,  , drop = FALSE]

  # build genmap_list in the SAME chromosome order used above
  chr_levels <- unique(map2$chr)  # already ordered by the 'order(chr,pos)'
  genmap_list <- lapply(chr_levels, function(cc) {
    as.matrix(map2$pos[map2$chr == cc])
  })

  # One-trait shortcut: covariance=TRUE doesn't add anything
  if (ncol(U) == 1L) covariance <- FALSE

  # ---- Threads (match parameter name n.Threads) ----
  if (length(n.Threads) != 1L || !is.finite(n.Threads) || n.Threads < 1 || n.Threads != as.integer(n.Threads)) {
    stop("`n.Threads` must be a positive integer.")
  }
  nThreads <- as.integer(n.Threads)

  # ---- Call C++ (Allier) ----
  temp <- cpp_calculate_covariance_RIL_allier(
    Crosses   = crosses,
    genMap    = genmap_list,
    M         = M,
    U         = U,
    t         = as.numeric(t),
    intensity = intensity,
    covariance = covariance,
    calcgains = calculate.gains,
    gains = gains,
    nThreads  = nThreads
  )

  # ---- Format outputs ----
  name_vec <- paste0(rep(c("EGBV","var","SPV"), each = ncol(U)), seq_len(ncol(U)))

  if (covariance) {
    if(calculate.gains){
      temp1 <- as.data.frame(temp$cross_values)[,1:(3*ncol(U))]
      names(temp1) <- name_vec
      temp2 <- as.data.frame(temp$cross_values)[,((3*ncol(U))+1):ncol(as.data.frame(temp$cross_values))]
      names(temp2) <- c("IDG_A","VARIDG_A","SPVIDG_A")
      return(list(cross_values=temp1,gains=temp2,covariances=temp$covariances))
    }else{
      temp1 <- as.data.frame(temp$cross_values)[1:(3*ncol(U))]
      names(temp1) <- name_vec
      return(list(cross_values=temp1,covariances=temp$covariances))
    }

  } else {
    out <- as.data.frame(temp)
    names(out) <- name_vec

    return(out)
  }
}
