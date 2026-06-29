test_that("poet_block preserves identity matrices", {
  out <- poet_block(diag(8), regularizer = "linear_shrinkage")

  expect_lt(max(abs(out$R_hat - diag(8))), 1e-8)
  expect_equal(length(out$clusters), 8)
  expect_equal(out$n_clusters, 8)
})

test_that("poet_block supports all residual regularizers", {
  set.seed(2)
  X <- matrix(rnorm(80 * 12), 80, 12)
  R <- cor(X)
  regularizers <- c("thresholding", "tapering", "linear_shrinkage", "nonlinear_shrinkage")

  for (regularizer in regularizers) {
    out <- if (identical(regularizer, "nonlinear_shrinkage")) {
      poet_block(R, X = X, n = 80, regularizer = regularizer)
    } else {
      poet_block(R, n = 80, regularizer = regularizer)
    }
    values <- eigen(out$R_hat, symmetric = TRUE, only.values = TRUE)$values

    expect_equal(dim(out$R_hat), c(12L, 12L))
    expect_true(isTRUE(all.equal(out$R_hat, t(out$R_hat), tolerance = 1e-8)))
    expect_lt(max(abs(diag(out$R_hat) - 1)), 1e-8)
    expect_gt(min(values), -1e-8)
    expect_equal(out$regularizer, regularizer)
  }
})

test_that("poet_block returns cluster metadata", {
  R <- diag(12)
  R[1:4, 1:4] <- 0.9
  R[5:8, 5:8] <- 0.85
  R[9:12, 9:12] <- 0.8
  diag(R) <- 1

  out <- poet_block(R, K_min = 3)

  expect_equal(length(out$clusters), 12)
  expect_equal(length(out$r_per_cluster), out$n_clusters)
  expect_gte(out$n_clusters, 3)
})
