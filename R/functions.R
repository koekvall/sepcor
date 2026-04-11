#' Estimate the covariance matrix Sigma = D (C2 \%x\% C1) D using the multivariate
#' normal likelihood, where \%x\% denotes the Kronecker product.
#'
#' @param E Matrix of dimension (n_rows n_cols) x n_obs, each column is a residual vector
#' @param n_rows Number of "rows" of the inverse vectorized columns of E
#' @param sepcov If TRUE, assume separable covariance matrix
#' @param tol Algorithm terminates when the relative increase in log-likelihood is less than tol
#' @param maxiter The maximum number of iterations the algorithm runs if not converging before
#' @param verbose Print additional info about iterates if TRUE
#' @param lambda Nuclear norm penalty parameter (>= 0). Subtracts (lambda/2)*tr(Sigma^{-1})
#'   from the log-likelihood, shrinking C1 and C2 toward the identity (independence).
#'   Only used when sepcov = FALSE. Default 0 gives the MLE.
#' @param n_starts Number of random starting points. The first start uses the default
#'   initialization (C1 = I, C2 = I). Additional starts use random correlation matrices
#'   and perturbed standard deviations. Only used when sepcov = FALSE. Default 1.
#' @return Final iterates, log-likelihood evaluated at these iterates, iterations,
#' and convergence info
#'
#' @export
#' @useDynLib sepcor
#' @importFrom Rcpp evalCpp
#' @importFrom stats rWishart
sepcor <- function(E, n_rows, sepcov = FALSE, tol = 1e-16, maxiter = 1000,
  verbose = FALSE, lambda = 0, n_starts = 1L)
{
  if(!is.matrix(E)){stop("E needs to be a rc x n matrix of residuals")}
  n_obs <- ncol(E)
  if(!(is.atomic(n_rows) && length(n_rows) == 1)){stop("n_rows needs to be a positive integer")}
  if(n_rows  <=0 || (n_rows != floor(n_rows))){stop("n_rows needs to be a positive integer")}
  n_cols <- nrow(E) / n_rows
  if(n_cols != floor(n_cols)){stop("Incompatible dimensions of E and n_rows")}
  if(!is.logical(sepcov)){stop("sepcov needs to TRUE or FALSE")}
  if(!(is.atomic(tol) && length(tol) == 1)){stop("tol needs to be a positive scalar")}
  if(tol < 0){stop("tol needs to be a positive scalar")}
  if(!(is.atomic(maxiter))){stop("maxiter needs to be a positive integer")}
  if(!(length(maxiter) == 1 && maxiter > 0 && (maxiter == floor(maxiter)))){stop("maxiter needs to be a positive integer")}
  if(!is.logical(verbose)){stop("verbose needs to be TRUE or FALSE")}
  if(lambda < 0){stop("lambda must be non-negative")}
  n_starts <- as.integer(n_starts)
  if(n_starts < 1L){stop("n_starts must be a positive integer")}

  if(sepcov){
    fit <- sepcov_rcpp(E, n_rows, tol, maxiter, verbose)
    return(fit)
  }

  # Helper: generate a random correlation matrix by rescaling a Wishart draw
  rand_corr <- function(d){
    W <- rWishart(1, d + 1, diag(d))[,,1]
    D_inv <- diag(1 / sqrt(diag(W)))
    D_inv %*% W %*% D_inv
  }

  S <- tcrossprod(E) / n_obs
  q <- nrow(E)

  # First start: default initialization
  best_fit <- sepcor_rcpp(E, diag(S), n_rows, tol, maxiter, verbose, lambda)

  # Additional random starts
  if(n_starts > 1L){
    for(s in 2:n_starts){
      # Random W: perturb sample standard deviations
      W_init <- sqrt(diag(S)) * exp(rnorm(q, 0, 0.5))
      fit_s <- sepcor_rcpp(E, W_init^2, n_rows, tol, maxiter, FALSE, lambda)
      if(fit_s$ll > best_fit$ll){
        best_fit <- fit_s
      }
    }
  }
  return(best_fit)
}

#' Compute standard errors for a fitted separable correlation model via
#' the expected Fisher information.
#'
#' @param fit A list returned by \code{sepcor} (with sepcov = FALSE).
#' @param E Matrix of dimension (n_rows * n_cols) x n_obs of residual vectors.
#' @param n_rows Number of rows of the inverse vectorized columns of E.
#' @return A list with components:
#'   \item{se_C2}{Standard errors for the upper-triangular entries of C2.}
#'   \item{se_C1}{Standard errors for the upper-triangular entries of C1.}
#'   \item{se_D}{Standard errors for the diagonal entries of D.}
#'   \item{vcov}{The full asymptotic covariance matrix of all parameters.}
#'
#' @export
sepcor_se <- function(fit, E, n_rows)
{
  n_obs <- ncol(E)
  nr    <- n_rows
  nc    <- nrow(E) / nr
  q     <- nr * nc

  C2  <- fit$C2
  C1  <- fit$C1
  C2i <- solve(C2)
  C1i <- solve(C1)

  # Upper-triangle indices for C2 and C1, matching upper.tri() column-major order
  U_ut <- which(upper.tri(C2), arr.ind = TRUE)   # n_U x 2
  V_ut <- which(upper.tri(C1), arr.ind = TRUE)   # n_V x 2
  n_U  <- nrow(U_ut)
  n_V  <- nrow(V_ut)

  aU <- U_ut[, 1L]; bU <- U_ut[, 2L]   # row/col indices for C2 upper tri
  aV <- V_ut[, 1L]; bV <- V_ut[, 2L]   # row/col indices for C1 upper tri

  C2i_ut <- C2i[U_ut]   # C2i entries at upper-tri positions
  C1i_ut <- C1i[V_ut]

  # For each diagonal element of D (index m, 1-indexed, column-major):
  #   m1 = row index in C1, m2 = column index in C2
  m1 <- ((seq_len(q) - 1L) %% nr) + 1L
  m2 <- ((seq_len(q) - 1L) %/% nr) + 1L

  # Block (C2, C2): I[j,k] = n*nr*(C2i[aj,ak]*C2i[bj,bk] + C2i[aj,bk]*C2i[bj,ak])
  I_UU <- n_obs * nr * (C2i[aU, aU] * C2i[bU, bU] +
                         C2i[aU, bU] * C2i[bU, aU])

  # Block (C1, C1): I[j,k] = n*nc*(C1i[aj,ak]*C1i[bj,bk] + C1i[aj,bk]*C1i[bj,ak])
  I_VV <- n_obs * nc * (C1i[aV, aV] * C1i[bV, bV] +
                         C1i[aV, bV] * C1i[bV, aV])

  # Block (C2, C1): I[j,k] = 2*n*C2i[aj,bj]*C1i[ak,bk]
  I_UV <- 2 * n_obs * outer(C2i_ut, C1i_ut)

  di <- 1 / as.vector(fit$D)   # reciprocal diagonal entries of D

  # Block (C2, D): I[j,m] = n*di[m]*C2i[aj,bj] * (1[m2==aj] + 1[m2==bj])
  I_UD <- n_obs * outer(C2i_ut, di) * (outer(aU, m2, "==") + outer(bU, m2, "=="))

  # Block (C1, D): I[j,m] = n*di[m]*C1i[aj,bj] * (1[m1==aj] + 1[m1==bj])
  I_VD <- n_obs * outer(C1i_ut, di) * (outer(aV, m1, "==") + outer(bV, m1, "=="))

  # Block (D, D): I[k,l] = n*(di[k]^2*delta[k,l] + di[k]*di[l]*(C2*C2i)[k2,l2]*(C1*C1i)[k1,l1])
  I_DD <- n_obs * (diag(di^2) + outer(di, di) * kronecker(C2 * C2i, C1 * C1i))

  # Assemble full information matrix in parameter order [C2 upper tri, C1 upper tri, D]
  I_full <- rbind(
    cbind(I_UU,    I_UV,    I_UD),
    cbind(t(I_UV), I_VV,    I_VD),
    cbind(t(I_UD), t(I_VD), I_DD)
  )

  vcov <- tryCatch(solve(I_full), error = function(e) {
    warning("Information matrix is singular; standard errors may be unreliable.")
    matrix(NA_real_, nrow(I_full), ncol(I_full))
  })

  se <- sqrt(pmax(diag(vcov), 0))

  list(
    se_C2   = se[seq_len(n_U)],
    se_C1   = se[n_U + seq_len(n_V)],
    se_D    = se[n_U + n_V + seq_len(q)],
    vcov    = vcov
  )
}

prof_log_lik <- function(Sigma_c, E)
#' Calculate (profile) log-likelihood for Sigma_c, i.e. the multivariate normal likelihood 
#' proportional to: -log(determinant(Sigma)) - trace{t(E) solve(Sigma) E},
#'  where Sigma =  Sigma_c t(Sigma_c).
#' 
#' @param Sigma_c An rc x rc lower triangular Cholesky factor to evaluate the likelihood at
#' @param E An rc x n_obs matrix of residual vectors, where n_obs is the number of such vectors
#' @return The likelihood value evaluated at the arguments supplied, *including constants*
#' 
#' @export
#' @useDynLib sepcor
#' @importFrom Rcpp evalCpp
{
  if(!(is.matrix(Sigma_c) && ncol(Sigma_c) == nrow(Sigma_c) && sum(Sigma_c[upper.tri(Sigma_c)]) == 0 &&
    sum(diag(Sigma_c) < 0)) == 0){stop("Sigma_c should be a lower triangular Cholesky root")}
  if(!(is.matrix(E) && nrow(E) == nrow(Sigma_c))){stop("E should be an rc x n matrix of residuals")}
  return(prof_log_lik_rcpp(Sigma_c, E))
}