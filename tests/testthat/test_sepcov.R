n_obs <- 1e5
n_rows <- 3
n_cols <- 3
n_resp <- n_rows * n_cols
U0 <- 0.5^abs(outer(1:n_cols, 1:n_cols, "-"))
V0 <- 0.5^abs(outer(1:n_rows, 1:n_rows, "-"))


Sig0_c <- kronecker(t(chol(U0)), t(chol(V0)))
Sig0 <- tcrossprod(Sig0_c)

Y <- Sig0_c %*% matrix(rnorm(n_resp * n_obs), ncol = n_obs)

fit <- sepcor(Y, n_rows, sepcov = TRUE)

testthat::test_that("algorithm converged", {
  testthat::expect_equal(fit$info, 0)
})

Sig_hat <- kronecker(fit$U, fit$V)

testthat::test_that("estimates seem consistent", {
  testthat::expect_lt(max(abs(Sig0 - Sig_hat)), 0.01)
})
