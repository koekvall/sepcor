#ifndef CONVENIENCE_H
#define CONVENIENCE_H

#include <RcppArmadillo.h>

arma::mat chol_solve(arma::mat&, arma::mat);

arma::mat tri_solve(arma::mat&, arma::mat, char, char);

extern "C"{ 
	void dpotrs_(char* UPLU, int* N, int* NRHS, double* A, int* LDA,
	double* B, int* LDB, int* INFO);
/* UPLU is 'U' if upper tri is stored, 'L' if lower tri
N is order of A
NRHS is number of columns of B
A is triangular Cholesky factor of A
LDA is leading dimension of A
B is the right hand side matrix on entry, solution on exit
LDB is leading dimension of B
INFO is 0 if success and -i if ith argument had illegal value*/
}
#endif