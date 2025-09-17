#' Tapering Regularization for LD Matrix (with configurable k_gap)
#'
#' Apply tapering regularization to a symmetric matrix.
#' For a given bandwidth k, entries with |i - j| <= k are kept as 1,
#' entries with k < |i - j| <= k_gap * k are linearly tapered to 0,
#' and entries with |i - j| > k_gap * k are set to 0.
#' The function searches over k_vec and returns the first matrix whose
#' minimum eigenvalue is at least eigenmin (else it returns the last one).
#'
#' @param A        A symmetric matrix (e.g., LD matrix).
#' @param k_vec    A vector of candidate bandwidths (default: 2:10).
#' @param k_gap    Taper reach multiplier (>= 1). k_gap = 2 means taper to 0 by 2k; k_gap = 3 means by 3k. Default: 2.
#' @param eigenmin Minimum eigenvalue threshold (default: 0.01).
#' @param print    Logical; if TRUE, prints chosen k and min eigenvalue.
#' @importFrom CppMatrix matrixEigen
#' @return         The regularized matrix after tapering.
#' @export
tapering <- function(A, k_vec = seq(2, 10), k_gap = 2, eigenmin = 0.01, print = FALSE) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  stopifnot(is.numeric(k_gap), length(k_gap) == 1L, k_gap >= 1)
  p <- nrow(A)

  # Build Toeplitz weight matrix given k and k_gap
  make_weight <- function(k, p, k_gap) {
    d <- 0:(p - 1)
    w <- numeric(length(d))

    # 1) Inside the band: weight = 1
    w[d <= k] <- 1

    # 2) Taper region: k < d <= k_gap * k
    upper_raw <- k_gap * k
    upper <- min(ceiling(upper_raw), p - 1)

    if (upper > k) {
      taper_idx <- which(d > k & d <= upper)
      # Linear decay from 1 at d = k to 0 at d = k_gap * k
      # Handle k_gap == 1 edge case (no taper region)
      denom <- (upper_raw - k)
      if (denom > 0) {
        w[taper_idx] <- (upper_raw - d[taper_idx]) / denom
      } else {
        # k_gap == 1 -> no taper; already handled by band
        w[taper_idx] <- 0
      }
    }

    # 3) Outside: 0 (already initialized)
    toeplitz(w)
  }

  # Helper: min eigenvalue with CppMatrix fallback to base eigen
  min_eig <- function(M) {
    val <- tryCatch(CppMatrix::matrixEigen(M)$values, error = function(e) NULL)
    if (is.null(val)) {
      val <- tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
    }
    if (all(is.na(val))) return(NA_real_)
    val[length(val)]
  }

  for (k in k_vec) {
    W <- make_weight(k, p, k_gap)
    Areg <- A * W
    Areg <- (Areg + t(Areg)) / 2

    lam_min <- min_eig(Areg)

    if (!is.na(lam_min) && lam_min >= eigenmin) {
      if (print) message("Selected k = ", k, " with min eigenvalue = ", round(lam_min, 6))
      return(Areg)
    }
  }

  # If none met the threshold, use the last k
  k_last <- tail(k_vec, 1L)
  W <- make_weight(k_last, p, k_gap)
  Areg <- A * W
  Areg <- (Areg + t(Areg)) / 2
  lam_min <- min_eig(Areg)
  if (print) message("No k met eigenmin; using k = ", k_last,
                     " with min eigenvalue = ", round(lam_min, 6))
  Areg
}
