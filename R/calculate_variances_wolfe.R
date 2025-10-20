#' Additive and dominance segregation variance of families from specified crosses (Wolfe)
#'
#' Calculate expected genomic estimated breeding value (EGBV), expected total genetic value (ETGV),
#' additive and dominance segregation variance, superior progeny value with additive variance (SPV),
#' and superior progeny value including dominance variance (TSPV) for F1 families.
#'
#' @param crosses A data.frame or matrix (n_crosses x 2) of parental indices
#'   (row indices into \code{hap1}/\code{hap2}) specifying the crosses.
#' @param hap1 Numeric haplotype matrix (individuals x markers) for the first haplotype,
#'   entries in \code{c(0, 1)} for \code{c("A", "B")}.
#' @param hap2 Numeric haplotype matrix (individuals x markers) for the second haplotype,
#'   entries in \code{c(0, 1)} for \code{c("A", "B")}.
#' @param genetic.map A data.frame with columns:
#'   \itemize{
#'     \item \code{site}: integer marker index (1..ncol(hap1))
#'     \item \code{chr}: chromosome identifier (numeric or factor)
#'     \item \code{pos}: position on the chromosome in cM (numeric)
#'   }
#'   All markers in \code{hap1}/\code{hap2} must appear in \code{genetic.map$site}.
#' @param U Numeric matrix of additive marker effects (markers in rows, traits in columns).
#'   Must have \code{nrow(U) == ncol(hap1)}.
#' @param D Numeric matrix of dominance marker effects (markers in rows, traits in columns).
#'   Must have \code{nrow(D) == ncol(hap1)} and \code{ncol(D) == ncol(U)}.
#' @param intensity Double. Standardized selection differential (used for SPV).
#' @param covariance Logical. If \code{TRUE}, also compute the segregation
#'   covariance matrix between all supplied traits.
#' @param calculate.gains Logical. If \code{TRUE}, calculate the desired gains index
#'  based on segregation (co)variances. This only works if covariance is TRUE, else will be ignored.
#' @param gains a vector of length equal to the number of traits with values representing the desired gains
#'
#' @param nThreads Integer (default 4). Number of OpenMP threads (if enabled at compile time).
#'
#' @return If \code{covariance = TRUE}, a list with element \code{cross_values}
#'   (a data.frame) whose columns are named \code{EGBV1, ETGV1, varA1, SPV1, varD, SPTV1,  EGBV2, ...} and a list with element \code{covariances} with segregation covariance matrix for each cross.
#'   If \code{covariance = TRUE}, a list with element \code{cross_values}
#'   (a data.frame) whose columns are named \code{EGBV1, ETGV1, varA1, SPV1, varD, SPTV1,  EGBV2, ...}, element \code{gains} (a data.frame) with \code{EGBVDG, ETGVDG, varADG, SPVDG, varDDG, SPTVDG} for the desired gains index
#'   and a list with element \code{covariances} with segregation covariance matrix for each cross.
#'   If \code{covariance = FALSE}, a data.frame with the same columns \code{EGBV1, ETGV1, varA1, SPV1, varD, SPTV1,  EGBV2, ....}.
#'
#' @examples
#' \dontrun{
#' # toy shapes only
#' crosses <- matrix(c(1,2, 3,4), ncol = 2, byrow = TRUE)
#' hap1 <- matrix(sample(0:1, 20, TRUE), nrow = 5)  # 5 inds x 4 markers
#' hap2 <- matrix(sample(0:1, 20, TRUE), nrow = 5)
#' genetic.map <- data.frame(site = 1:ncol(hap1), chr = c(1,1,2,2), pos = c(0,10,0,10))
#' U <- matrix(rnorm(ncol(hap1)), ncol = 1)
#' D <- matrix(rnorm(ncol(hap1)), ncol = 1)
#' out <- calculate_variances_F1(crosses, genetic.map, hap1, hap2, U, D, intensity = 1.0, covariance = FALSE)
#' }
#' @export
calculate_variances_F1 <- function(crosses, genetic.map, hap1, hap2, U, D, intensity, covariance,calculate.gains=FALSE, gains=NULL, nThreads = 4L) {

  hap1 <- as.matrix(hap1); hap2 <- as.matrix(hap2)
  U <- as.matrix(U); D <- as.matrix(D)
  crosses <- as.matrix(crosses)

  #  Basic shape checks
  if (!is.logical(covariance) || length(covariance) != 1L) {
    stop("`covariance` must be a single logical (TRUE/FALSE).")
  }
  if (!is.numeric(intensity) || length(intensity) != 1L) {
    stop("`intensity` must be a single numeric (standardized selection differential).")
  }
  if (!is.numeric(nThreads) || length(nThreads) != 1L || nThreads < 1) {
    stop("`nThreads` must be a positive integer.")
  }

  if (ncol(hap1) <= 0L) stop("hap1 must have markers in columns.")
  if (!identical(dim(hap1), dim(hap2))) stop("hap1 and hap2 must have identical dimensions.")
  if (nrow(U) != ncol(hap1)) {
    stop("U must have nrow(U) == ncol(hap1). Found: nrow(U) = ", nrow(U), ", ncol(hap1) = ", ncol(hap1), ".")
  }
  if (nrow(D) != ncol(hap1)) {
    stop("D must have nrow(D) == ncol(hap1). Found: nrow(D) = ", nrow(D), ", ncol(hap1) = ", ncol(hap1), ".")
  }
  if (ncol(D) != ncol(U)) stop("U and D must have the same number of trait columns.")

  if (ncol(crosses) != 2L) stop("`crosses` must have exactly 2 columns (P1, P2).")

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

  #  Validate crosses against haplotypes
  genos <- unique(as.vector(crosses))
  if (any(!is.finite(genos))) stop("`crosses` contains non-finite entries.")
  if (any(genos < 1 | genos > nrow(hap1))) {
    stop("Some genotype indices in `crosses` are outside 1..nrow(hap1).")
  }

  #  Validate and align genetic.map with markers
  req_cols <- c("site", "chr", "pos")
  if (!all(req_cols %in% names(genetic.map))) {
    stop("`genetic.map` must contain columns: ", paste(req_cols, collapse = ", "), ".")
  }
  if (any(genetic.map$site < 1 | genetic.map$site > ncol(hap1))) {
    stop("`genetic.map$site` contains indices outside 1..ncol(hap1).")
  }
  if (!all(seq_len(ncol(hap1)) %in% genetic.map$site)) {
    stop("Some markers in hap1/hap2 are missing from `genetic.map$site`.")
  }

  # Reorder marker columns of haplotypes and effects to match the map order
  ord <- as.integer(genetic.map$site)
  hap1 <- hap1[, ord, drop = FALSE]
  hap2 <- hap2[, ord, drop = FALSE]
  U    <- U[ord, , drop = FALSE]
  D    <- D[ord, , drop = FALSE]

  #  Build per-chromosome position matrices for C++
  chrs <- unique(genetic.map$chr)
  genmap_list <- lapply(chrs, function(x) {
    as.matrix(genetic.map[genetic.map$chr == x, "pos", drop = TRUE])
  })

  #  Call C++
  temp <- cpp_calculate_covariance_wolfe(
    Crosses   = crosses,
    genMap    = genmap_list,
    Hap1      = hap1,
    Hap2      = hap2,
    U         = U,
    D         = D,
    intensity = intensity,
    covariance = covariance,
    calcgains = calculate.gains,
    gains = gains,
    nThreads  = as.integer(nThreads)
  )

  #  Format outputs
  name_vec <- paste0(rep(c("EGBV","ETGV","var_a","SPV","var_d","TSPV"), each = ncol(U)), seq_len(ncol(U)))

  if (covariance) {
    if (calculate.gains) {
      temp1 <- as.data.frame(temp$cross_values)[ , 1:(6*ncol(U)), drop = FALSE]
      names(temp1) <- name_vec
      temp2 <- as.data.frame(temp$cross_values)[ , (6*ncol(U) + 1):(6*ncol(U) + 6), drop = FALSE]
      names(temp2) <- c("IDG_A","VARIDG_A","SPVIDG_A","IDG_AD","VARIDG_AD","SPVIDG_AD")
      return(list(cross_values = temp1, gains = temp2, additive_covariances = temp$covA, dominance_covariances = temp$covD))
    } else {
      temp1 <- as.data.frame(temp$cross_values)[ , 1:(6*ncol(U)), drop = FALSE]
      names(temp1) <- name_vec
      return(list(cross_values = temp1, additive_covariances = temp$covA, dominance_covariances = temp$covD))
    }
  } else {
    out <- as.data.frame(temp)
    names(out) <- name_vec
    return(out)
  }

}
