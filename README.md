# LDRegularization

The `LDRegularization` package provides a collection of covariance and correlation matrix regularization methods, with a focus on high-dimensional settings such as linkage disequilibrium (LD) matrices in statistical genetics. The package implements sparse regularization, shrinkage, nonlinear shrinkage, POET-style estimators, and covariance transfer learning helpers.

These tools help ensure positive semi-definiteness and improve estimation accuracy when working with large, noisy covariance or correlation matrices.

## Installation

You can install the `LDRegularization` package from GitHub using the `devtools` package:

```r
devtools::install_github("harryyiheyang/LDRegularization")
```

## Functions

The package currently includes the following main functions:

- poet_thresholding: POET estimator with entrywise thresholding on the idiosyncratic covariance.

- poet_banding: POET estimator with banding (keep entries within a fixed bandwidth).

- poet_tapering: POET estimator with tapering (apply a tapering kernel depending on |i-j|).

- poet_linear_shrinkage: POET estimator with linear shrinkage of the idiosyncratic covariance.

- poet_nonlinear_shrinkage: POET estimator with mixed nonlinear shrinkage of
  the idiosyncratic covariance.

- thresholding: Standalone entrywise MCP thresholding for covariance/correlation matrices.

- banding: Standalone banding with bandwidth `K`.

- tapering: Standalone tapering with bandwidth `K`.

- linear_shrinkage: Standalone linear shrinkage.

- nonlinear_shrinkage: Standalone mixed nonlinear shrinkage for individual-level
  data.

- covariance_tl: Covariance transfer learning helper with covariance or
  projection source transfer.

- select_tl_alpha: 5-fold cross-validation helper for choosing the transfer
  strength `tl_alpha`.

Most matrix estimators accept either a precomputed matrix (`S` or `A`) or
individual-level data `X`. Nonlinear shrinkage requires `X`. The POET functions
automatically select the number of latent factors via the ratio-type criterion.
Sparse methods use theory-rate defaults only when the user does not supply the
tuning parameter. Standalone thresholding uses
`lambda = 2 * sqrt(log(p) / n)`, while POET thresholding uses
`lambda = 2 * max(sqrt(log(p) / n), 1 / sqrt(p))`. Banding and tapering use
`K = ceiling(n^(1 / (2 * alpha + 1)))` with `alpha = 1` by default, capped at
`floor(p / 2)`. These defaults require `n`; otherwise pass `lambda` or `K`
directly. Standalone linear shrinkage uses `alpha = 0.05` by default; POET
linear shrinkage uses `alpha = 0.5` by default. Nonlinear shrinkage has
`shrinkage = 0` by default, so scripts should set it explicitly, for example
`shrinkage = 0.5`, when a stronger nonlinear component is desired.

Sparse estimators are made positive definite with fixed-support positive-definite (FSPD) linear shrinkage, using `eigenmin = 0.001` by default. FSPD is a final safeguard and does not choose the sparsity tuning parameter.

Covariance transfer learning uses target correlation `S` and source
information only to guide the factor space. With `strategy = "covariance"`, the
factor space is estimated from `S + tl_alpha * source`. With
`strategy = "projection"`, the source enters through a trace-normalized source
basis/projector. The residual component is then regularized by one existing
method: thresholding, banding, tapering, linear shrinkage, or mixed nonlinear
shrinkage. `select_tl_alpha` chooses `tl_alpha` from a small grid using target
individual-level data `X_target` and off-diagonal Frobenius loss on held-out
target correlations. If `X_source` is supplied and no source covariance is
given, the source correlation is computed internally.

## Dependencies

The package makes use of efficient matrix operations implemented in CppMatrix, which relies on Rcpp and RcppArmadillo for performance.

## Example

```r
library(LDRegularization)

# Example: POET with thresholding
set.seed(123)
p <- 50
n <- 200
X <- matrix(rnorm(n * p), n, p)
S <- cov(X)
Sigma_hat <- poet_thresholding(S, n = n)
```

## License

This package is licensed under the MIT License.

## Contact

Yihe Yang
Email: yxy1234@case.edu
