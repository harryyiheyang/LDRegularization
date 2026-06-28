#' Covariance transfer learning regularization
#'
#' Use source covariance information to guide the factor space, then regularize
#' the target idiosyncratic component with one of the existing sparse methods.
#'
#' @param S Optional target covariance or correlation matrix.
#' @param source Optional source information. For `strategy = "covariance"`, this must
#'   be a covariance or correlation matrix. For `strategy = "projection"`, this
#'   can be a source basis matrix with `p` rows, a source projector, or a source
#'   covariance matrix from which a basis is extracted.
#' @param n Target sample size used to form `S`.
#' @param X Optional target individual-level data matrix. Used to compute `S`
#'   when `S` is not supplied and required for
#'   `method = "nonlinear_shrinkage"`.
#' @param R_source Optional source covariance or correlation matrix.
#' @param X_source Optional source individual-level data. If `source` and
#'   `R_source` are omitted, `R_source = matrixCor(X_source)` is used.
#' @param P_source Optional source basis or projection matrix for
#'   `strategy = "projection"`.
#' @param strategy Source transfer strategy: `"covariance"` or `"projection"`.
#'   `"P"` is accepted as an alias for `"projection"`.
#' @param method Residual regularization method.
#' @param tl_alpha Nonnegative transfer strength. `0` means no transfer.
#' @param lambda Optional threshold for `method = "thresholding"`.
#' @param K Optional bandwidth for `method = "banding"` or `"tapering"`.
#' @param alpha Smoothness parameter used only for default `K`.
#' @param shrink_alpha Shrinkage intensity used only for
#'   `method = "linear_shrinkage"`.
#' @param nonlinear_shrinkage Mixing weight used only for
#'   `method = "nonlinear_shrinkage"`. Default is 0.
#' @param eigenmin Minimum eigenvalue target for FSPD. Default is 0.001.
#' @param source_rank Optional rank used when `strategy = "projection"` and
#'   `source` is a square covariance matrix.
#'
#' @return A positive-definite transfer-learning regularized correlation
#'   matrix.
#' @export
covariance_tl <- function(
    S = NULL,
    source = NULL,
    n = NULL,
    X = NULL,
    R_source = NULL,
    X_source = NULL,
    P_source = NULL,
    strategy = c("covariance", "projection", "P"),
    method = c(
      "thresholding",
      "banding",
      "tapering",
      "linear_shrinkage",
      "nonlinear_shrinkage"
    ),
    tl_alpha = 0.5,
    lambda = NULL,
    K = NULL,
    alpha = 1,
    shrink_alpha = 0.5,
    nonlinear_shrinkage = 0,
    eigenmin = 1e-3,
    source_rank = NULL
) {
  strategy <- .ld_match_tl_strategy(strategy)
  method <- match.arg(method)
  input <- .ld_resolve_input(S = S, X = X, n = n, name = "S")
  S <- .ld_clean_correlation(input$S, "S")
  source_info <- .ld_resolve_tl_source_input(
    source = source,
    R_source = R_source,
    X_source = X_source,
    P_source = P_source,
    strategy = strategy,
    p = nrow(S),
    source_rank = source_rank
  )

  .ld_fit_covariance_tl_prepared(
    S = S,
    X = input$X,
    n = input$n,
    source_info = source_info,
    strategy = strategy,
    tl_alpha = tl_alpha,
    method = method,
    lambda = lambda,
    K = K,
    alpha = alpha,
    shrink_alpha = shrink_alpha,
    nonlinear_shrinkage = nonlinear_shrinkage,
    eigenmin = eigenmin
  )
}

#' Select transfer strength for covariance transfer learning
#'
#' Select `tl_alpha` by K-fold cross-validation on target individual-level data.
#' The validation score is the off-diagonal Frobenius loss on correlation scale.
#'
#' @param X_target Target individual-level data matrix with observations in rows.
#' @param source Source information passed to `covariance_tl()`.
#' @param R_source Optional source covariance or correlation matrix.
#' @param X_source Optional source individual-level data. If `source` and
#'   `R_source` are omitted, `R_source = matrixCor(X_source)` is used.
#' @param P_source Optional source basis or projection matrix for
#'   `strategy = "projection"`.
#' @param strategy Source transfer strategy: `"covariance"` or `"projection"`.
#'   `"P"` is accepted as an alias for `"projection"`.
#' @param method Residual regularization method.
#' @param tl_alpha_grid Candidate nonnegative transfer strengths. `0` is added
#'   automatically if not supplied.
#' @param folds Number of folds. Default is 5.
#' @param lambda Optional threshold for `method = "thresholding"`.
#' @param K Optional bandwidth for `method = "banding"` or `"tapering"`.
#' @param alpha Smoothness parameter used only for default `K`.
#' @param shrink_alpha Shrinkage intensity used only for
#'   `method = "linear_shrinkage"`.
#' @param nonlinear_shrinkage Mixing weight used only for
#'   `method = "nonlinear_shrinkage"`. Default is 0.
#' @param eigenmin Minimum eigenvalue target for FSPD. Default is 0.001.
#' @param source_rank Optional rank used when `strategy = "projection"` and
#'   `source` is a square covariance matrix.
#'
#' @return A list with the selected `tl_alpha`, average CV scores, fold scores,
#'   and the fold assignments.
#' @export
select_tl_alpha <- function(
    X_target,
    source = NULL,
    R_source = NULL,
    X_source = NULL,
    P_source = NULL,
    strategy = c("covariance", "projection", "P"),
    method = c(
      "thresholding",
      "banding",
      "tapering",
      "linear_shrinkage",
      "nonlinear_shrinkage"
    ),
    tl_alpha_grid = c(0, 0.05, 0.1, 0.25, 0.5, 1, 2),
    folds = 5,
    lambda = NULL,
    K = NULL,
    alpha = 1,
    shrink_alpha = 0.5,
    nonlinear_shrinkage = 0,
    eigenmin = 1e-3,
    source_rank = NULL
) {
  if (!is.matrix(X_target)) {
    X_target <- as.matrix(X_target)
  }
  storage.mode(X_target) <- "double"
  n_total <- nrow(X_target)
  p <- ncol(X_target)
  if (n_total < 4L) {
    stop("X_target must have at least 4 rows for cross-validation.", call. = FALSE)
  }

  folds <- as.integer(round(folds[1L]))
  if (!is.finite(folds) || folds < 2L) {
    stop("folds must be at least 2.", call. = FALSE)
  }
  folds <- min(folds, floor(n_total / 2))

  tl_alpha_grid <- as.numeric(tl_alpha_grid)
  tl_alpha_grid <- tl_alpha_grid[is.finite(tl_alpha_grid) & tl_alpha_grid >= 0]
  if (!length(tl_alpha_grid)) {
    stop("tl_alpha_grid must contain at least one nonnegative finite value.", call. = FALSE)
  }
  tl_alpha_grid <- sort(unique(c(0, tl_alpha_grid)))

  strategy <- .ld_match_tl_strategy(strategy)
  method <- match.arg(method)
  source_info <- .ld_resolve_tl_source_input(
    source = source,
    R_source = R_source,
    X_source = X_source,
    P_source = P_source,
    strategy = strategy,
    p = p,
    source_rank = source_rank
  )

  fold_id <- sample(rep(seq_len(folds), length.out = n_total))
  fold_scores <- matrix(
    NA_real_,
    nrow = folds,
    ncol = length(tl_alpha_grid),
    dimnames = list(paste0("fold", seq_len(folds)), paste0("alpha_", tl_alpha_grid))
  )
  fold_sizes <- integer(folds)

  for (fold in seq_len(folds)) {
    val_idx <- which(fold_id == fold)
    train_idx <- which(fold_id != fold)
    fold_sizes[fold] <- length(val_idx)

    R_train <- .ld_matrix_cor(X_target[train_idx, , drop = FALSE])
    R_val <- .ld_matrix_cor(X_target[val_idx, , drop = FALSE])

    for (j in seq_along(tl_alpha_grid)) {
      estimate <- .ld_fit_covariance_tl_prepared(
        S = R_train,
        X = X_target[train_idx, , drop = FALSE],
        n = length(train_idx),
        source_info = source_info,
        strategy = strategy,
        tl_alpha = tl_alpha_grid[j],
        method = method,
        lambda = lambda,
        K = K,
        alpha = alpha,
        shrink_alpha = shrink_alpha,
        nonlinear_shrinkage = nonlinear_shrinkage,
        eigenmin = eigenmin
      )
      fold_scores[fold, j] <- .ld_offdiag_frobenius_loss(estimate, R_val)
    }
  }

  cv_score <- as.numeric(crossprod(fold_sizes, fold_scores) / sum(fold_sizes))
  best_idx <- which.min(cv_score)
  scores <- data.frame(
    tl_alpha = tl_alpha_grid,
    score = cv_score
  )

  list(
    tl_alpha = tl_alpha_grid[best_idx],
    scores = scores,
    fold_scores = fold_scores,
    folds = fold_id,
    strategy = strategy,
    method = method
  )
}
