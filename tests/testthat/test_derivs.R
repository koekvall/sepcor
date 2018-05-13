n_obs <- 30
n_rows <- 3
n_cols <- 3
n_resp <- n_rows * n_cols
U0 <- 0.5^abs(outer(1:n_cols, 1:n_cols, "-"))
V0 <- 0.5^abs(outer(1:n_rows, 1:n_rows, "-"))


Sig0_c <- kronecker(t(chol(U0)), t(chol(V0)))
Sig0 <- tcrossprod(Sig0_c)

Y <- Sig0_c %*% matrix(rnorm(n_resp * n_obs), ncol = n_obs)

# Testing U and V gradient
log_lik <- function(theta){
  U <- matrix(0, nrow = 3, ncol = 3)
  V <- matrix(0, 3, 3)
  U[upper.tri(U, diag = TRUE)] <- theta[1:6]
  U[lower.tri(U)] <- t(U)[lower.tri(U)]
  V[upper.tri(V, diag = TRUE)] <- theta[7:12]
  V[lower.tri(V)] <- t(V)[lower.tri(V)]
  Sig_c <- kronecker(t(chol(U)), t(chol(V)))
  return(sepcor::prof_log_lik(E = Y, Sigma_c = Sig_c))
}

theta0 <- c(as.vector(U0[upper.tri(U0, diag = TRUE)]), as.vector(V0[upper.tri(V0, diag = TRUE)]))
numerical_grad <- numDeriv::grad(log_lik, theta0)

gradUV <- sepcor:::ll_grad_UV(Y, U0, V0)
U_grad <- 2 * matrix(gradUV[1:9], 3, 3)
diag(U_grad) <- diag(U_grad) / 2
V_grad <- 2 * matrix(gradUV[10:18], 3, 3)
diag(V_grad) <- diag(V_grad) / 2
analytic_grad <- c(U_grad[upper.tri(U_grad, diag = T)], 
  V_grad[upper.tri(V_grad, diag = T)]) 

testthat::test_that("gradient for U, V is correct", {
  testthat::expect_equal(analytic_grad, numerical_grad)
})

# Testing W gradient
log_lik <- function(theta){
  Sig_c <- sweep(kronecker(t(chol(U0)), t(chol(V0))), 1, 1 / theta, "*")
  return(sepcor::prof_log_lik(E = Y, Sigma_c = Sig_c))
}

numerical_grad <- numDeriv::grad(log_lik, rep(1, 9))
analytic_grad <- as.vector(sepcor:::ll_grad_W_inv(W = rep(1, 9), E = Y, U = U0, V = V0))
testthat::test_that("gradient for W is correct", {
  testthat::expect_equal(analytic_grad, numerical_grad)
})
