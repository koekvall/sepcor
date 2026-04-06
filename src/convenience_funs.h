#ifndef CONVENIENCE_H
#define CONVENIENCE_H

#include <RcppArmadillo.h>

arma::mat chol_solve(arma::mat&, arma::mat);

arma::mat tri_solve(arma::mat&, arma::mat, char, char);

#endif
