#include <RcppArmadillo.h>
#include "convenience_funs.h"

arma::mat chol_solve(arma::mat& A, arma::mat B)
{
	// argument A is the lower Cholesky root of p.d. matrix in the system (A A^T)X = B
	// Solve L X1 = B, then L^T X = X1
	arma::mat X = arma::solve(arma::trimatl(A), B, arma::solve_opts::fast);
	return arma::solve(arma::trimatu(A.t()), X, arma::solve_opts::fast);
}

arma::mat tri_solve(arma::mat& A, arma::mat B, char UPLO, char TRANS)
{
  // argument A is triangular A in the system AX = B (or A^T X = B if TRANS = 'T')
  if(TRANS == 'T' || TRANS == 't'){
    if(UPLO == 'L' || UPLO == 'l'){
      return arma::solve(arma::trimatu(A.t()), B, arma::solve_opts::fast);
    } else {
      return arma::solve(arma::trimatl(A.t()), B, arma::solve_opts::fast);
    }
  } else {
    if(UPLO == 'L' || UPLO == 'l'){
      return arma::solve(arma::trimatl(A), B, arma::solve_opts::fast);
    } else {
      return arma::solve(arma::trimatu(A), B, arma::solve_opts::fast);
    }
  }
}