#include <RcppArmadillo.h>
#include "convenience_funs.h"

arma::mat chol_solve(arma::mat& A, arma::mat B)
{
	// argument A is the lower Cholesky root of p.d. matrix A in the system AX = B,
  // stored as a regular matrix.
	char UPLO = 'L';
  int N = A.n_rows;
	int NRHS = B.n_cols;
	int LDB = B.n_rows;
	int INFO;
	dpotrs_(&UPLO, &N, &NRHS, A.memptr(), &N, B.memptr(), &LDB, &INFO);

  if(INFO != 0){
    Rcpp::stop("Cholesky solve failed");
  }
	return B;
}

arma::mat tri_solve(arma::mat& A, arma::mat B, char UPLO, char TRANS)
{
  // argument A is triangular A in the system AX = B
  char DIAG = 'N';
  int N = A.n_rows;
  int NRHS = B.n_cols;
  int LDB = B.n_rows;
  int INFO;
  arma::dtrtrs_(&UPLO, &TRANS, &DIAG, &N, &NRHS, A.memptr(), &N, B.memptr(), &LDB, &INFO);

  if(INFO != 0){
    Rcpp::stop("Triangular solve failed");
  }
  return B;
}