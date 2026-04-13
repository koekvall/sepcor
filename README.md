# sepcor

R package for likelihood-based inference with separable correlation matrices.
The covariance matrix is modeled as

Sigma = D (C2 ⊗ C1) D,

where C1 and C2 are correlation matrices, D is a diagonal matrix of positive
standard deviations, and ⊗ denotes the Kronecker product. Unlike separable
covariance, the variances in D are unrestricted.

## Installation

```r
# install.packages("devtools")
devtools::install_github("koekvall/sepcor")
```

## Usage

```r
library(sepcor)

# E is a (r * c) x n matrix of residual vectors
fit <- sepcor(E, n_rows = r)
fit$C1   # estimated row correlation matrix (r x r)
fit$C2   # estimated column correlation matrix (c x c)
fit$D    # estimated standard deviations (r*c vector)
fit$ll   # log-likelihood at the estimate

# Standard errors from expected Fisher information
se <- sepcor_se(fit, E, n_rows = r)
se$se_C1  # SEs for off-diagonal entries of C1
se$se_C2  # SEs for off-diagonal entries of C2
se$se_D   # SEs for diagonal entries of D
```

## Functions

- `sepcor()` — Fits the model by maximum likelihood via block-coordinate
  ascent. Supports regularization (`lambda`) and multiple random starts
  (`n_starts`). Set `sepcov = TRUE` to fit the separable covariance model.
- `sepcor_se()` — Computes standard errors from the expected Fisher
  information.
- `prof_log_lik_sep()` — Evaluates the profile log-likelihood for given
  C1, C2, D.
- `prof_log_lik()` — Evaluates the profile log-likelihood for an
  unstructured covariance (via its Cholesky factor).

## Reference

Ekvall, K. O. (2025+). Likelihood-based inference with separable correlation
matrices. *Stat*, to appear.

## Author

Karl Oskar Ekvall (k.ekvall@ufl.edu)

## License

MIT
