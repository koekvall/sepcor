#include <RcppArmadillo.h>
#include "convenience_funs.h"
#include "likelihood_funs.h"

//[[Rcpp::export]]
double prof_log_lik_rcpp(arma::mat& Sigma_chol, const arma::mat& E)
{
	const unsigned int q = Sigma_chol.n_rows;
	const unsigned int n = E.n_cols;

	double ll = -0.5 * n * q * std::log(2 * arma::datum::pi);
	ll -=  n * arma::accu(arma::log(Sigma_chol.diag()));
	ll -= 0.5 * arma::accu(arma::square(tri_solve(Sigma_chol, E, 'L', 'N')));
	return ll;
}


//[[Rcpp::export]]
arma::vec ll_grad_UV(arma::mat& E, arma::mat U, arma::mat V)
{
  const unsigned int n_rows = V.n_cols;
  const unsigned int n_cols = U.n_cols; 
  const unsigned int n_resp = E.n_rows;
  const unsigned int n_obs = E.n_cols;

  arma::mat U_grad = (n_obs * n_rows) * U; //n_obs because of scaling after loop
  arma::mat V_grad = (n_obs * n_cols) * V;
  U = arma::chol(U, "lower");
  V = arma::chol(V, "lower");
  for(size_t ii = 0; ii < n_obs; ii++){
    arma::mat Ei(E.colptr(ii), n_rows, n_cols, false, true);
    U_grad -= Ei.t() * chol_solve(V, Ei);
    V_grad -= Ei * chol_solve(U, Ei.t());
  }
  // Pre and post multiply by U^{-1}
  U_grad = chol_solve(U, U_grad);
  U_grad = chol_solve(U, U_grad.t()).t();
  // Pre and post multiply by V^{-1}
  V_grad = chol_solve(V, V_grad);
  V_grad = chol_solve(V, V_grad.t()).t();
  return -0.5 * arma::join_cols(arma::vectorise(U_grad), arma::vectorise(V_grad));
}

//[[Rcpp::export]]
arma::vec ll_grad_W_inv(arma::vec& W, arma::mat& E, arma::mat U, arma::mat V)
{
  arma::mat R_c = arma::kron(arma::chol(U, "lower"), arma::chol(V, "lower"));
  arma::mat W_inv_grad = chol_solve(R_c, arma::diagmat(arma::pow(W, -1.0)) *
    E * E.t());
  return -W_inv_grad.diag() + E.n_cols * W;
}