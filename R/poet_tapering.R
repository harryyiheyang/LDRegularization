#' POET with tapering and FSPD positive-definite modification
#'
#' Estimate a factor component from the leading eigenpairs of `S`, taper the
#' idiosyncratic covariance, and apply FSPD linear shrinkage to the sparse
#' idiosyncratic estimator.
#'
#' @param S Symmetric covariance or correlation matrix.
#' @param n Sample size used to form `S`.
#' @param K Optional bandwidth. Overrides the theoretical default.
#' @param alpha Bandable covariance smoothness parameter used in the theoretical
#'   bandwidth. Default is 1.
#' @param eigenmin Minimum eigenvalue target for FSPD. Default is 0.001.
#'
#' @return A positive-definite POET regularized correlation matrix.
#' @export
poet_tapering <- function(
    S,
    n,
    K = NULL,
    alpha = 1,
    eigenmin = 1e-3
) {
  n <- .ld_validate_n(n)
  comp <- .ld_poet_components(
    S,
    n,
    cutoff_method = "D.ratio",
    k_min = 5,
    k_max = min(15, floor(nrow(S) / 2))
  )

  E_reg <- tapering(
    comp$E,
    n = n,
    K = K,
    alpha = alpha,
    eigenmin = eigenmin
  )

  Shat <- .ld_symmetrize(comp$P + E_reg)
  stats::cov2cor(Shat)
}
