#' Estimate the covariance matrix Sigma = W k(U, V) W using the multivariate
#' normal likelihood, where k(U, V) means the Kronecker product of U and V.
#'
#' @param E Matrix of dimension (n_rows n_cols) x n_obs, each column is a residual vector
#' @param n_rows Number of "rows" of the inverse vectorized columns of E
#' @param sepcov If TRUE, assume separable covariance matrix
#' @param tol Algorithm terminates when an iteration increases the log-likelihood less than tol
#' @param maxiter The maximum number of iterations the algorithm runs if not converging before
#' @param verbose Print additional info about iterates if TRUE
#' @return Final iterates, log-likelihood evaluated at these iterates, iterations,
#' and convergence info 
#'
#' @export
#' @useDynLib sepcor
#' @importFrom Rcpp evalCpp
sepcor <- function(E, n_rows, sepcov = FALSE, tol = 1e-16, maxiter = 1000, 
  verbose = FALSE)
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

  if(!sepcov){
    S <- tcrossprod(E) / n_obs
    fit <- sepcor_rcpp(E, diag(S), n_rows, tol, maxiter, verbose)
  } else{
    fit <- sepcov_rcpp(E, n_rows, tol, maxiter, verbose)
  }
  return(fit)
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