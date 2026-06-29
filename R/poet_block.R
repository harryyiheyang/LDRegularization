#' Block-wise POET for LD matrix regularization
#'
#' Regularize an LD correlation matrix by clustering variants on absolute
#' correlation, building a block-wise POET factor space, refining that space
#' with a few subspace iterations, and regularizing the residual with one of
#' four simple methods.
#'
#' @param R Symmetric LD correlation or covariance matrix.
#' @param n Optional sample size. Used for theory-rate defaults when available.
#' @param regularizer Residual regularization method. One of
#'   `"thresholding"`, `"tapering"`, `"linear_shrinkage"`, or
#'   `"nonlinear_shrinkage"`.
#' @param X Optional individual-level data matrix. Required when
#'   `regularizer = "nonlinear_shrinkage"`.
#' @param r_threshold Minimum pairwise absolute correlation used by the
#'   fallback clustering cut. Default is 0.95.
#' @param cumvar Per-cluster cumulative variance proportion for eigenvector
#'   retention. Default is 0.95.
#' @param K_min Minimum number of clusters. Capped at `ncol(R)`.
#' @param n_iter Number of subspace iteration steps. Default is 3.
#' @param linkage Hierarchical clustering linkage method. Default is
#'   `"complete"`.
#' @param eigenmin Minimum eigenvalue target for positive-definite correction.
#'   Default is 0.001.
#' @param lambda Optional threshold for `regularizer = "thresholding"`.
#' @param K Optional bandwidth for `regularizer = "tapering"`.
#' @param alpha Optional shrinkage intensity for
#'   `regularizer = "linear_shrinkage"`, or smoothness parameter for the
#'   default bandwidth in tapering.
#'
#' @return A list with components `R_hat`, `H`, `Lambda`, `E_hat`, `clusters`,
#'   `r_per_cluster`, `n_clusters`, `regularizer`, and `tuning`.
#' @export
poet_block <- function(
    R,
    n = NULL,
    regularizer = c("thresholding", "tapering", "linear_shrinkage", "nonlinear_shrinkage"),
    X = NULL,
    r_threshold = 0.95,
    cumvar = 0.95,
    K_min = 5,
    n_iter = 3,
    linkage = "complete",
    eigenmin = 1e-3,
    lambda = NULL,
    K = NULL,
    alpha = NULL
) {
  regularizer <- match.arg(regularizer)
  linkage <- match.arg(
    linkage,
    c("complete", "single", "average", "mcquitty", "median", "centroid", "ward.D", "ward.D2")
  )

  R <- .ld_as_square_matrix(R, "R")
  if (nrow(R) < 2L) {
    stop("R must have at least two rows and columns.", call. = FALSE)
  }
  if (anyNA(R)) {
    stop("R must not contain missing values.", call. = FALSE)
  }
  R <- .ld_clean_correlation(R, "R")

  p <- nrow(R)
  n <- if (is.null(n)) NULL else .ld_validate_n(n)
  r_threshold <- .ld_scalar_in_range(r_threshold, "r_threshold", lower = 0, upper = 1)
  cumvar <- .ld_scalar_in_range(cumvar, "cumvar", lower = 0, upper = 1)
  K_min <- min(max(as.integer(round(K_min[1L])), 1L), p)
  n_iter <- max(as.integer(round(n_iter[1L])), 0L)
  if (length(eigenmin) != 1L || !is.finite(eigenmin) || eigenmin <= 0) {
    stop("eigenmin must be a positive finite scalar.", call. = FALSE)
  }

  clusters <- .ld_poet_block_clusters(
    R,
    r_threshold = r_threshold,
    K_min = K_min,
    linkage = linkage
  )
  basis <- .ld_poet_block_basis(R, clusters = clusters, cumvar = cumvar)
  factor <- .ld_poet_block_subspace(R, basis$U, n_iter = n_iter)

  E <- .ld_poet_block_residual_covariance(R - factor$R_low)
  residual <- .ld_poet_block_regularize_residual(
    E,
    regularizer = regularizer,
    X = X,
    H = factor$H,
    n = n,
    lambda = lambda,
    K = K,
    alpha = alpha,
    eigenmin = eigenmin
  )

  R_hat <- .ld_symmetrize(factor$R_low + residual$E_hat)
  R_hat <- .ld_fspd(R_hat, eigenmin = eigenmin)
  R_hat <- .ld_clean_correlation(R_hat, "R_hat")
  R_hat <- .ld_fspd(R_hat, eigenmin = eigenmin)
  R_hat <- .ld_clean_correlation(R_hat, "R_hat")

  list(
    R_hat = R_hat,
    H = factor$H,
    Lambda = factor$Lambda,
    E_hat = residual$E_hat,
    clusters = clusters,
    r_per_cluster = basis$r_per_cluster,
    n_clusters = max(clusters),
    regularizer = regularizer,
    tuning = residual$tuning
  )
}

.ld_scalar_in_range <- function(x, name, lower, upper) {
  if (length(x) != 1L || !is.finite(x) || x <= lower || x > upper) {
    stop(name, " must be in (", lower, ", ", upper, "].", call. = FALSE)
  }
  as.numeric(x)
}

.ld_poet_block_clusters <- function(R, r_threshold, K_min, linkage) {
  p <- nrow(R)
  D <- 1 - abs(R)
  D[!is.finite(D)] <- 1
  D[D < 0] <- 0
  diag(D) <- 0

  hc <- stats::hclust(stats::as.dist(D), method = linkage)
  heights <- sort(hc$height)
  clusters <- NULL

  if (length(heights) > 1L) {
    gaps <- diff(heights)
    candidates <- order(gaps, decreasing = TRUE)
    for (idx in candidates) {
      h_cut <- heights[idx] + 1e-10
      candidate <- stats::cutree(hc, h = h_cut)
      if (length(unique(candidate)) >= K_min) {
        clusters <- candidate
        break
      }
    }
  }

  if (is.null(clusters)) {
    clusters <- stats::cutree(hc, h = 1 - r_threshold)
  }
  if (length(unique(clusters)) < K_min) {
    clusters <- stats::cutree(hc, k = K_min)
  }

  as.integer(factor(clusters, levels = sort(unique(clusters))))
}

.ld_poet_block_basis <- function(R, clusters, cumvar) {
  p <- nrow(R)
  n_clusters <- max(clusters)
  vectors <- vector("list", n_clusters)
  r_per_cluster <- integer(n_clusters)

  for (cluster_id in seq_len(n_clusters)) {
    rows <- which(clusters == cluster_id)
    R_cluster <- R[rows, rows, drop = FALSE]
    eig <- .ld_eigen(R_cluster)
    values <- pmax(eig$values, 0)
    total <- sum(values)
    if (!is.finite(total) || total <= .Machine$double.eps) {
      r_cluster <- 1L
    } else {
      r_cluster <- which(cumsum(values) / total >= cumvar)[1L]
      r_cluster <- min(max(as.integer(r_cluster), 1L), length(rows))
    }
    vectors[[cluster_id]] <- eig$vectors[, seq_len(r_cluster), drop = FALSE]
    r_per_cluster[cluster_id] <- r_cluster
  }

  total_rank <- sum(r_per_cluster)
  U <- matrix(0, p, total_rank)
  col_start <- 1L
  for (cluster_id in seq_len(n_clusters)) {
    rows <- which(clusters == cluster_id)
    cols <- col_start:(col_start + r_per_cluster[cluster_id] - 1L)
    U[rows, cols] <- vectors[[cluster_id]]
    col_start <- max(cols) + 1L
  }

  list(U = U, r_per_cluster = r_per_cluster)
}

.ld_poet_block_subspace <- function(R, H, n_iter) {
  if (n_iter > 0L) {
    for (i in seq_len(n_iter)) {
      H <- .ld_matrix_multiply(R, H)
      H <- qr.Q(qr(H))[, seq_len(ncol(H)), drop = FALSE]
    }
  }

  RH <- .ld_matrix_multiply(R, H)
  Lambda <- .ld_symmetrize(crossprod(H, RH))
  R_low <- .ld_matrix_multiply(.ld_matrix_multiply(H, Lambda), t(H))
  list(H = H, Lambda = Lambda, R_low = .ld_symmetrize(R_low))
}

.ld_poet_block_regularize_residual <- function(
    E,
    regularizer,
    X,
    H,
    n,
    lambda,
    K,
    alpha,
    eigenmin
) {
  p <- nrow(E)

  if (identical(regularizer, "linear_shrinkage")) {
    shrinkage <- if (is.null(alpha)) 0.5 else alpha
    if (length(shrinkage) != 1L || !is.finite(shrinkage)) {
      stop("alpha must be a finite scalar.", call. = FALSE)
    }
    shrinkage <- min(max(shrinkage[1L], 0), 1)
    target <- diag(diag(E), p, p)
    E_hat <- .ld_symmetrize((1 - shrinkage) * E + shrinkage * target)
    return(list(E_hat = E_hat, tuning = list(alpha = shrinkage)))
  }

  if (identical(regularizer, "thresholding")) {
    resolved_lambda <- if (is.null(lambda)) {
      .ld_poet_block_default_threshold(p, n)
    } else {
      lambda
    }
    E_hat <- thresholding(E, lambda = resolved_lambda, eigenmin = eigenmin)
    return(list(E_hat = E_hat, tuning = list(lambda = resolved_lambda)))
  }

  if (identical(regularizer, "nonlinear_shrinkage")) {
    if (is.null(X)) {
      stop(
        "X is required when regularizer = \"nonlinear_shrinkage\".",
        call. = FALSE
      )
    }
    E_hat <- .ld_nonlinear_residual(
      X,
      H,
      shrinkage = 0,
      eigenmin = eigenmin,
      scale = FALSE
    )
    return(list(E_hat = E_hat, tuning = list(shrinkage = 0)))
  }

  smoothness <- if (is.null(alpha)) 1 else alpha
  resolved_K <- if (is.null(K)) {
    .ld_poet_block_default_bandwidth(p, n, smoothness)
  } else {
    K
  }

  E_hat <- tapering(E, K = resolved_K, alpha = smoothness, eigenmin = eigenmin)

  list(E_hat = E_hat, tuning = list(K = resolved_K, alpha = smoothness))
}

.ld_poet_block_residual_covariance <- function(E) {
  E <- .ld_fix_residual_diag(E)
  d <- diag(E)
  bad <- !is.finite(d) | d <= .Machine$double.eps
  if (any(bad)) {
    positive_min <- suppressWarnings(min(d[!bad], na.rm = TRUE))
    if (!is.finite(positive_min)) {
      positive_min <- .Machine$double.eps
    }
    d[bad] <- positive_min
    diag(E) <- d
  }
  .ld_symmetrize(E)
}

.ld_poet_block_default_threshold <- function(p, n) {
  if (!is.null(n)) {
    return(.ld_theory_threshold(p, n))
  }
  2 * sqrt(log(max(p, 2L)) / max(p, 2L))
}

.ld_poet_block_default_bandwidth <- function(p, n, alpha) {
  if (!is.null(n)) {
    return(.ld_theory_bandwidth(p, n, alpha))
  }
  .ld_theory_bandwidth(p, max(p, 2L), alpha)
}

.ld_matrix_multiply <- function(x, y) {
  tryCatch(
    CppMatrix::matrixMultiply(x, y),
    error = function(e) x %*% y
  )
}
