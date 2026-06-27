#' Adaptive linear shrinkage for a correlation matrix
#'
#' Shrink a sample correlation matrix toward a target matrix. If `lambda` is not
#' supplied, a Ledoit-Wolf-style plug-in value is estimated from the centered
#' data matrix `X`.
#'
#' @param S Symmetric correlation matrix.
#' @param X Centered data matrix with the same number of columns as `S`.
#' @param lambda Optional shrinkage intensity in `[0, 1]`.
#' @param target Optional shrinkage target. Defaults to the identity matrix.
#'
#' @return A regularized correlation matrix after shrinkage.
#' @export
linear_shrinkage <- function(S, X, lambda = NULL, target = NULL) {
  S <- .ld_as_square_matrix(S, "S")
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (ncol(X) != ncol(S)) {
    stop("X must have the same number of columns as S.", call. = FALSE)
  }

  p <- ncol(S)
  if (is.null(target)) {
    target <- diag(1, p)
  } else {
    target <- .ld_as_square_matrix(target, "target")
    if (!all(dim(target) == c(p, p))) {
      stop("target must have the same dimensions as S.", call. = FALSE)
    }
  }

  S <- .ld_symmetrize(S)
  target <- .ld_symmetrize(target)

  if (is.null(lambda)) {
    row_norm4 <- rowSums(X^2)^2
    pi_hat <- mean(row_norm4) - sum(S^2)
    d2_hat <- sum((S - target)^2)
    lambda <- if (d2_hat > 0) pi_hat / d2_hat else 1
  }

  lambda <- min(max(lambda[1L], 0), 1)

  .ld_symmetrize(lambda * target + (1 - lambda) * S)
}
