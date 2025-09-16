#' Adaptive Linear Shrinkage for Correlation Matrix
#'
#' Apply linear shrinkage to a correlation matrix.
#' If lambda is not provided, estimate it adaptively from data (Ledoit–Wolf).
#'
#' @param S A symmetric correlation matrix with diagonal elements equal to 1.
#' @param X Data matrix (n × p), assumed centered (mean zero by column).
#' @param lambda Shrinkage intensity in [0,1]. If NULL, estimate adaptively.
#' @param target Target matrix for shrinkage (default: identity).
#' @param print Logical; if TRUE, print chosen lambda.
#' @return A regularized correlation matrix after shrinkage.
#' @export
linear_shrinkage <- function(S,X,lambda = NULL,target = NULL,print = FALSE) {
stopifnot(is.matrix(S), nrow(S) == ncol(S))
stopifnot(is.matrix(X), ncol(X) == ncol(S))

n <- nrow(X)
p <- ncol(X)
if (is.null(target)) target <- diag(1, p)

## symmetrize
S <- (S + t(S)) / 2

## estimate lambda if NULL
if (is.null(lambda)) {
row_norm4 <- rowSums(X^2)^2
pi_hat <- mean(row_norm4) - sum(S^2)

d2_hat <- sum((S - target)^2)

lambda <- max(0, min(1, pi_hat / d2_hat))
if (print) {
message("Estimated lambda = ", round(lambda, 4))
}
}

## shrinkage
Sigma_hat <- lambda * target + (1 - lambda) * S
return(Sigma_hat)
}
