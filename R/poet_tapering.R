#' poet_tapering: POET with tapering on the idiosyncratic covariance
#'
#' Factor+Sparse (POET) regularization where the factor part uses the top
#' eigenpairs of S and the idiosyncratic part is tapered entrywise.
#'
#' @param S               A symmetric matrix (p×p): covariance OR correlation.
#' @param n               Sample size used to estimate S.
#' @param cutoff_method   "D.ratio" (default) or "ratio" for picking #factors m.
#' @param k_min,k_max     Search range for m. Defaults: 5 and min(15, floor(p/2)).
#' @param tapering_vec    Candidate taper bandwidths (passed to `tapering()`).
#' @param k_gap           Taper reach multiplier in `tapering()` (2 ⇒ to 2k, 3 ⇒ to 3k).
#' @param eigenmin        Minimum eigenvalue floor used inside `tapering()`. Default 1e-3.
#' @param return_correlation Logical; if TRUE return correlation, else covariance. Default TRUE.
#'
#' @return A PSD matrix (p×p): correlation if return_correlation=TRUE, else covariance.
#' @importFrom CppMatrix matrixEigen matrixMultiply
#' @export
poet_tapering <- function(
    S, n,
    cutoff_method = "D.ratio",
    k_min = 5,
    k_max = min(15, floor(nrow(S) / 2)),
    tapering_vec = seq(2, 20, 2),
    k_gap = 2,
    eigenmin = 1e-3,
    return_correlation = TRUE
) {
  stopifnot(is.matrix(S), nrow(S) == ncol(S))
  p <- ncol(S)

  # Symmetrize (no hard-thresholding to avoid breaking PSD)
  S <- (S + t(S)) / 2

  # Eigendecomposition (descending)
  eig <- tryCatch(CppMatrix::matrixEigen(S), error = function(e) eigen(S, symmetric = TRUE))
  d <- as.numeric(eig$values)
  U <- eig$vectors
  # Safety: ensure decreasing order
  ord <- order(d, decreasing = TRUE)
  d <- d[ord]; U <- U[, ord, drop = FALSE]

  # Bound search range
  k_min <- max(2L, as.integer(k_min))
  k_max <- min(as.integer(k_max), p - 2L)
  if (k_min > k_max) k_min <- max(2L, min(k_max, p - 2L))

  # Select number of factors m (pck)
  z <- rep(NA_real_, p)
  if (cutoff_method == "D.ratio") {
    j_lo <- max(k_min, 2L); j_hi <- min(k_max, p - 1L)
    if (j_lo <= j_hi) {
      for (j in j_lo:j_hi) {
        den <- (d[j] - d[j + 1L])
        num <- (d[j - 1L] - d[j])
        z[j - 1L] <- if (den != 0) num / den else NA_real_
      }
    }
  } else if (cutoff_method == "ratio") {
    j_lo <- max(k_min, 2L); j_hi <- min(k_max, p)
    if (j_lo <= j_hi) {
      for (j in j_lo:j_hi) {
        z[j - 1L] <- if (d[j] != 0) d[j - 1L] / d[j] else NA_real_
      }
    }
  } else {
    stop("Unknown cutoff_method.")
  }
  idx <- suppressWarnings(which.max(z))
  pck <- if (length(idx) == 0L || is.na(idx)) max(k_min, 1L) else idx + 1L
  pck <- min(max(pck, 1L), p - 1L)

  # Factor component P = U_k diag(d_k) U_k^T  (use stable construction)
  Uk  <- U[, seq_len(pck), drop = FALSE]
  dk  <- d[seq_len(pck)]
  # Optional MCP shrink on factor eigenvalues (requires user-defined mcp())
  # Estimate idiosyncratic variance scale (POET bias correction style)
  trS <- sum(diag(S))
  denom <- (p - pck) - (pck * p / n)
  hatc <- if (denom > 0) (trS - sum(dk)) / denom else 0
  if (exists("mcp", mode = "function")) {
    dk <- mcp(dk, hatc * p / n, a = 3)
  }
  # Build P via scaling columns of U_k by sqrt(dk)
  Ukd <- Uk %*% diag(pmax(dk, 0), nrow = pck, ncol = pck)
  P   <- Uk %*% diag(pmax(dk, 0), pck, pck) %*% t(Uk)
  P   <- (P + t(P)) / 2

  # Idiosyncratic part
  E <- S - P
  E <- (E + t(E)) / 2
  # Ensure non-negative variances
  ee <- diag(E)
  pos_min <- suppressWarnings(min(ee[ee > 0]))
  if (is.finite(pos_min)) {
    ee[ee < 0] <- pos_min
  } else {
    ee[ee < 0] <- 0
  }
  diag(E) <- ee

  # Taper the idiosyncratic covariance
  E_tap <- tapering(E, k_vec = tapering_vec, k_gap = k_gap, eigenmin = eigenmin)

  # Recompose
  Shat <- P + E_tap
  Shat <- (Shat + t(Shat)) / 2

  if (isTRUE(return_correlation)) {
    return(cov2cor(Shat))
  } else {
    return(Shat)
  }
}
