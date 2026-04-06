#' Estimate the covariance matrix Sigma = W k(U, V) W using the multivariate
#' normal likelihood, where k(U, V) means the Kronecker product of U and V.
#'
#' @param E Matrix of dimension (n_rows n_cols) x n_obs, each column is a residual vector
#' @param n_rows Number of "rows" of the inverse vectorized columns of E
#' @param sepcov If TRUE, assume separable covariance matrix
#' @param tol Algorithm terminates when an iteration increases the log-likelihood less than tol
#' @param maxiter The maximum number of iterations the algorithm runs if not converging before
#' @param verbose Print additional info about iterates if TRUE
#' @param lambda Ridge penalty parameter (>= 0). Adds (lambda/2)[tr(U^{-1}) + tr(V^{-1})]
#'   to the negative log-likelihood, shrinking U and V toward the identity (independence).
#'   Only used when sepcov = FALSE. Default 0 gives the MLE.
#' @param n_starts Number of random starting points. The first start uses the default
#'   initialization (U = I, V = I). Additional starts use random correlation matrices
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
#' the observed Fisher information (numerical Hessian of the log-likelihood).
#'
#' @param fit A list returned by \code{sepcor} (with sepcov = FALSE).
#' @param E Matrix of dimension (n_rows * n_cols) x n_obs of residual vectors.
#' @param n_rows Number of rows of the inverse vectorized columns of E.
#' @return A list with components:
#'   \item{se_U}{Standard errors for the upper-triangular entries of U (row-major).}
#'   \item{se_V}{Standard errors for the upper-triangular entries of V (row-major).}
#'   \item{se_W}{Standard errors for the diagonal entries of W (on the log scale).}
#'   \item{vcov}{The full asymptotic covariance matrix of all parameters.}
#'
#' @export
sepcor_se <- function(fit, E, n_rows)
{
  if(!requireNamespace("numDeriv", quietly = TRUE)){
    stop("Package 'numDeriv' is needed for standard error computation.")
  }
  n_obs <- ncol(E)
  n_cols <- nrow(E) / n_rows

  U <- fit$U
  V <- fit$V
  W_vec <- as.vector(fit$W) # standard deviations

  # Indices for free parameters: upper triangle of U, upper triangle of V, log(W)
  U_ut <- which(upper.tri(U), arr.ind = TRUE)
  V_ut <- which(upper.tri(V), arr.ind = TRUE)
  n_U <- nrow(U_ut)
  n_V <- nrow(V_ut)
  q <- n_rows * n_cols

  # Pack MLE parameters into a single vector
  theta_hat <- c(
    U[upper.tri(U)],      # upper-triangular of U (off-diagonal correlations)
    V[upper.tri(V)],      # upper-triangular of V
    log(W_vec)            # log standard deviations
  )

  # Negative log-likelihood as a function of the free parameter vector
  nll <- function(theta){
    # Unpack U
    U_cur <- diag(n_cols)
    U_cur[upper.tri(U_cur)] <- theta[1:n_U]
    U_cur[lower.tri(U_cur)] <- t(U_cur)[lower.tri(U_cur)]

    # Unpack V
    V_cur <- diag(n_rows)
    V_cur[upper.tri(V_cur)] <- theta[(n_U + 1):(n_U + n_V)]
    V_cur[lower.tri(V_cur)] <- t(V_cur)[lower.tri(V_cur)]

    # Unpack W
    W_cur <- exp(theta[(n_U + n_V + 1):(n_U + n_V + q)])

    # Build Sigma and evaluate
    R <- kronecker(U_cur, V_cur)
    Sigma <- diag(W_cur) %*% R %*% diag(W_cur)

    # Check positive definiteness
    Sigma_chol <- tryCatch(chol(Sigma), error = function(e) NULL)
    if(is.null(Sigma_chol)) return(1e20)

    Sigma_c <- t(Sigma_chol) # lower triangular
    -prof_log_lik(Sigma_c, E)
  }

  H <- numDeriv::hessian(nll, theta_hat)

  # Invert to get asymptotic covariance (divide by n for the MLE)
  vcov <- tryCatch(solve(H), error = function(e){
    warning("Hessian is singular; standard errors may be unreliable.")
    matrix(NA, nrow(H), ncol(H))
  })

  se <- sqrt(pmax(diag(vcov), 0))

  list(
    se_U = se[1:n_U],
    se_V = se[(n_U + 1):(n_U + n_V)],
    se_logW = se[(n_U + n_V + 1):(n_U + n_V + q)],
    vcov = vcov
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