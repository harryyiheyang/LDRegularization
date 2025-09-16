#' Entrywise MCP thresholding on correlation; PSD check & return covariance
#'
#' @param S        Symmetric matrix; covariance or correlation (diag>=0).
#' @param lam_vec  Candidate lambdas (increasing), e.g. seq(0.01, 0.10, 0.01)
#' @param eigenmin Minimum eigenvalue target on covariance scale (default 0.01)
#' @param print    Whether to print chosen lambda
#' @param mcp_a    MCP 'a' parameter (default 3)
#' @importFrom CppMatrix matrixEigen
#' @return         Regularized covariance matrix (PSD up to eigenmin)
#' @export
thresholding <- function(
S,
lam_vec  = seq(0.01, 0.10, by = 0.01),
eigenmin = 0.01,
print    = FALSE,
mcp_a    = 3
) {
stopifnot(is.matrix(S), nrow(S) == ncol(S))
p <- nrow(S)
S <- (S + t(S)) / 2
diag_vals <- sqrt(pmax(diag(S), 0))
S_cor <- cov2cor(S)
diag(S_cor) <- 0

apply_mcp <- function(M, lam) {
x <- as.vector(M)
x_new <- mcp(x, lam = lam, a = mcp_a)
M_new <- matrix(x_new, nrow = p, ncol = p)
M_new <- (M_new + t(M_new)) / 2
diag(M_new) <- 1
M_new
}

best_Mcor <- NULL
best_lam  <- tail(lam_vec, 1L)
best_min  <- -Inf

for (lam in lam_vec) {
Mcor  <- apply_mcp(S_cor, lam)
Sigma <- t(t(Mcor) * diag_vals) * diag_vals
eigmin <- tryCatch({
ev <- matrixEigen((Sigma + t(Sigma)) / 2)$values
min(ev, na.rm = TRUE)
}, error = function(e) NA_real_)

if (!is.na(eigmin) && eigmin >= eigenmin) {
if (print) message("Selected lambda = ",lam," with min eigenvalue = ", round(eigmin, 6))
best_Mcor <- Mcor
best_lam  <- lam
break
}
if (!is.na(eigmin) && eigmin > best_min) {
best_min  <- eigmin
best_lam  <- lam
best_Mcor <- Mcor
}
}

if (is.null(best_Mcor)) best_Mcor <- apply_mcp(S_cor, tail(lam_vec, 1L))

Sigma <- t(t(best_Mcor) * diag_vals) * diag_vals
Sigma <- (Sigma + t(Sigma)) / 2
diag(Sigma) <- diag_vals^2

if (print && best_min > -Inf) {
message("No lambda met eigenmin; using lambda = ", best_lam,
    " with best min eigenvalue ≈ ", round(best_min, 6))
}
Sigma
}
