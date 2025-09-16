#' Banding Regularization for LD Matrix
#'
#' Apply banding regularization to a symmetric matrix.
#' For a given bandwidth parameter k, entries within |i-j| <= k are kept and
#' entries with |i-j| > k are set to 0. The function searches over a sequence
#' of k values (k_vec) and chooses the smallest k for which the minimum eigenvalue
#' of the regularized matrix is at least eigenmin.
#'
#' @param A A symmetric matrix (e.g., LD matrix).
#' @param k_vec A vector of candidate banding bandwidths (default: 2:10).
#' @param eigenmin Minimum eigenvalue threshold (default: 0.01).
#' @param print Logical; if TRUE, prints the chosen k and the minimum eigenvalue.
#' @importFrom CppMatrix matrixEigen
#' @return The regularized matrix after banding.
#' @export
banding <- function(A, k_vec = seq(2, 10), eigenmin = 0.01, print = FALSE) {
stopifnot(is.matrix(A), nrow(A) == ncol(A))
p <- nrow(A)

make_weight <- function(k, p) {
d <- 0:(p - 1)
w <- as.numeric(d <= k)
toeplitz(w)
}

for (k in k_vec) {
W <- make_weight(k, p)
Areg <- A * W
Areg <- (Areg + t(Areg)) / 2

lam_min <- tryCatch(
matrixEigen(Areg)$values[p],
error = function(e) NA_real_
)

if (!is.na(lam_min) && lam_min >= eigenmin) {
if (print) {
print(paste("Selected k =", k, "with min eigenvalue =", round(lam_min, 6)))
}
return(Areg)
}
}

k_last <- tail(k_vec, 1L)
W <- make_weight(k_last, p)
Areg <- A * W
Areg <- (Areg + t(Areg)) / 2
lam_min <- tryCatch(
matrixEigen(Areg)$values[p],
error = function(e) NA_real_
)
if (print) {
print(paste("No k met eigenmin; using k =", k_last,"with min eigenvalue =", round(lam_min, 6)))
}
return(Areg)
}
