n_obs <- 1e5
n_rows <- 3
n_cols <- 3
n_resp <- n_rows * n_cols
U0 <- 0.5^abs(outer(1:n_cols, 1:n_cols, "-"))
V0 <- 0.5^abs(outer(1:n_rows, 1:n_rows, "-"))

W0 <- sqrt(seq(1:n_resp) / n_resp)

Sig0_c <- sweep(kronecker(t(chol(U0)), t(chol(V0))), 1, W0, "*")
Sig0 <- tcrossprod(Sig0_c)

Y <- Sig0_c %*% matrix(rnorm(n_resp * n_obs), ncol = n_obs)

fit <- sepcor(Y, n_rows)

testthat::test_that("algorithm converged", {
  testthat::expect_equal(fit$info, 0)
})

Sig_hat <- sweep(sweep(kronecker(fit$U, fit$V), 1, fit$W, "*"), 2, fit$W, "*")


testthat::test_that("estimates seem consistent", {
  testthat::expect_lt(max(abs(Sig0 - Sig_hat)), 0.01)
})
