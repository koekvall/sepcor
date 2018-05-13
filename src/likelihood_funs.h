#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <RcppArmadillo.h>

double prof_log_lik_rcpp(arma::mat&, const arma::mat&);
arma::vec ll_grad_sepcov(arma::mat&, arma::mat, arma::mat);
#endif