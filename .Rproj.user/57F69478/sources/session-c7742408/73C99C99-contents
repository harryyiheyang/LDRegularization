# LDRegularization

The `LDRegularization` package provides a collection of covariance and correlation matrix regularization methods, with a focus on high-dimensional settings such as linkage disequilibrium (LD) matrices in statistical genetics. The package implements factor-based shrinkage estimators such as POET (Principal Orthogonal complEment Thresholding) with different strategies for regularizing the idiosyncratic component (thresholding, banding, tapering, and linear shrinkage).

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

- poet_tapering: POET estimator with tapering (apply a tapering kernel depending on |i−j|).

- poet_linear_shrinkage: POET estimator with linear shrinkage of the idiosyncratic covariance.

- thresholding: Standalone entrywise MCP thresholding for covariance/correlation matrices.

All functions automatically select the number of latent factors via ratio-type criteria and search over candidate regularization parameters (e.g., thresholding levels, banding widths, tapering bandwidths) to ensure the resulting covariance matrix is positive semi-definite.

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
