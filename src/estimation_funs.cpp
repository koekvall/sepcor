#include <RcppArmadillo.h>
#include "likelihood_funs.h"
#include "convenience_funs.h"

//[[Rcpp::export]]
Rcpp::List sepcor_rcpp(const arma::mat E, arma::vec W, const int n_rows,
 const double tol, const int maxiter, const bool verbose)
{
	// \Sigma = W (U \otimes V) W
	const unsigned int n_obs = E.n_cols;
	const unsigned int n_cols = E.n_rows / n_rows;

	arma::mat U = arma::eye(n_cols, n_cols);
	arma::mat U_c = U; // Cholesky factor storage
	arma::mat V = arma::eye(n_rows, n_rows);
	arma::mat V_c = V; // Cholesky factor storage

	const arma::mat S = E * E.t() * (1.0 / n_obs); // Residual covariance matrix
	arma::mat M(n_rows * n_cols, n_rows * n_cols); // Re-usable storage

	arma::mat E_tilde(n_rows, n_obs * n_cols); // \tilde{E}_i = W^{-1}E_i, i.e.
	// scaled residual vectors matricized, concatenated by columns

	int info = -1; // Keep track of convergence
	double ll_old = -DBL_MAX; // Initialize likelihood
	unsigned int iter = 0;

	W = arma::pow(W, -1.0); // use reciprocal standard deviations (precision)
	while(info != 0 && info != 1){
		bool success;
		bool converged;
		double ll_new;

		// Update U
		E_tilde = arma::reshape(arma::diagmat(W) * E, n_rows, n_obs * n_cols);
		U.zeros();
		for (size_t ii = 0; ii < n_obs; ++ii)
		{
			U += E_tilde.cols(ii * n_cols, (ii + 1) * n_cols - 1).t() * 
			chol_solve(V_c, E_tilde.cols(ii * n_cols, (ii + 1) * n_cols - 1));
		}
		U *= (1.0 / (n_obs * n_rows));
		// Rescale U to corr-mat and W accordingly
		U_c.col(0) = 1.0 / arma::sqrt(U.diag()); // Store precision param. updates
		for (int ii = 0; ii < W.n_elem; ++ii)
		{
			W(ii) *= U_c(ii / n_rows, 0);
		}
		U = arma::diagmat(U_c.col(0)) * U * arma::diagmat(U_c.col(0));
		success = arma::chol(U_c, U, "lower");
		if(!success){
			Rcpp::warning("Iterate of U not PD.");
			info = 2;
			break;
		}

		// Update V
		E_tilde = arma::reshape(arma::diagmat(W) * E, n_rows, n_obs * n_cols);
		V.zeros();
		for (size_t ii = 0; ii < n_obs; ++ii)
		{
			V += E_tilde.cols(ii * n_cols, (ii + 1) * n_cols - 1) *
				chol_solve(U_c, E_tilde.cols(ii * n_cols, (ii + 1) * n_cols - 1).t());
		}
		V *= (1.0 / (n_obs * n_cols));
		// Rescale V to corr-mat and W accordingly
		V_c.col(0) = 1.0 / arma::sqrt(V.diag()); // Store precision param. updates
		for (int ii = 0; ii < W.n_elem; ++ii)
		{
			W(ii) *= V_c(ii % n_rows , 0);
		}
		V = arma::diagmat(V_c.col(0)) * V * arma::diagmat(V_c.col(0));
		success = arma::chol(V_c, V, "lower");
		if(!success){
			Rcpp::warning("Iterate of V not PD.");
			info = 3;
			break;
		}

		// Update precision parameters
		M = S % arma::kron(chol_solve(U_c, arma::eye(n_cols, n_cols)), 
			chol_solve(V_c, arma::eye(n_rows, n_rows)));
		for (size_t ii = 0; ii < W.n_elem; ++ii)
		{
			double a = -M(ii, ii);
			double b = -arma::as_scalar(arma::accu(M.col(ii) % W) - W(ii) * M(ii, ii));
			W(ii) = (-b - std::sqrt(std::pow(b, 2.0) - 4.0 * a)) / (2.0 * a);
		}
		// Check convergence and wrap up iteration
		iter ++;
		if(iter >= maxiter){
			Rcpp::warning("Maximum iterations reached");
			info = 1;
		}
		
		if(iter % 50 == 0){
		  Rcpp::checkUserInterrupt();
		}
		
		M = arma::diagmat(arma::pow(W, -1.0)) * arma::kron(U_c, V_c); // use storage
		ll_new = prof_log_lik_rcpp(M, E);
		converged = ((ll_new -  ll_old)  < tol);
		if(converged){
			info = 0;
		}
		ll_old = ll_new;

		if(verbose){
			Rcpp::Rcout << "Iteration: " << iter << std::endl;
			Rcpp::Rcout << "loglik: " << ll_new << std::endl;
			Rcpp::Rcout << "U: " << U << std::endl;
			Rcpp::Rcout << "V: " << V << std::endl;
		}
	}
	return Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("V") = V,
		Rcpp::Named("W") = arma::pow(W, -1.0), Rcpp::Named("ll") = 
			ll_old, Rcpp::Named("iter") = iter, Rcpp::Named("info") = info);
}

//[[Rcpp::export]]
Rcpp::List sepcov_rcpp(arma::mat E, const int n_rows,
 const double tol, const int maxiter, const bool verbose)
{
	// \Sigma = U \otimes V
	const unsigned int n_obs = E.n_cols;
	const unsigned int n_cols = E.n_rows / n_rows;

	arma::mat U = arma::eye(n_cols, n_cols);
	arma::mat U_c = U; // Cholesky factor storage
	arma::mat V = arma::eye(n_rows, n_rows);
	arma::mat V_c = V;

	const arma::mat S = E * E.t() * (1.0 / n_obs);
	arma::mat M(n_rows * n_cols, n_rows * n_cols);

	// E in vec_inv form, same memory
	arma::mat E_tilde(E.memptr(), n_rows, n_cols * n_obs, false, true); 

	double ll_old = -DBL_MAX;
	unsigned int iter = 0;
	int info = -1; // Keep track of convergence
	while(info != 0 && info != 1){
	  bool success;
	  bool converged;
	  double ll_new;

		// Update U
		U.zeros();
		for (size_t ii = 0; ii < n_obs; ++ii){
			U += E_tilde.cols(ii * n_cols, (ii + 1) * n_cols - 1).t() * 
			chol_solve(V_c, E_tilde.cols(ii * n_cols, (ii + 1) * n_cols - 1));
		}
		U *= (1.0 / (n_obs * n_rows));
		success = arma::chol(U_c, U, "lower");
		if(!success){
			Rcpp::warning("Iterate of U not PD.");
			info = 2;
			break;
		}

		// Update V
		V.zeros();
		for (size_t ii = 0; ii < n_obs; ++ii){
			V += E_tilde.cols(ii * n_cols, (ii + 1) * n_cols - 1) * 
				chol_solve(U_c, E_tilde.cols(ii * n_cols, (ii + 1) * n_cols - 1).t());
		}
		V *= (1.0 / (n_obs * n_cols));
		success = arma::chol(V_c, V, "lower");
		if(!success){
			Rcpp::warning("Iterate of V not PD.");
			info = 3;
			break;
		}

		// Impose identifiability
		V *= U(0, 0);
		V_c *= std::sqrt(U(0, 0));
		U_c *= std::sqrt(1.0 / U(0, 0));
		U *= 1.0 / U(0, 0);

		// Check convergence and wrap up iteration
		iter ++;
		if(iter >= maxiter){
			Rcpp::warning("Maximum iterations reached");
			info = 1;
		}
		if(iter % 50 == 0){
		  Rcpp::checkUserInterrupt();
		}

		M = arma::kron(U_c, V_c);
		ll_new = prof_log_lik_rcpp(M, E);
		converged = ((ll_new -  ll_old)  < tol);
		if(converged){
			info = 0;
		}
		ll_old = ll_new;

		if(verbose){
			Rcpp::Rcout << "Iteration: " << iter << std::endl;
			Rcpp::Rcout << "loglik: " << ll_new << std::endl;
			Rcpp::Rcout << "U: " << U << std::endl;
			Rcpp::Rcout << "V: " << V << std::endl;
		}
	}
	return Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("V") = V,
		Rcpp::Named("ll") = ll_old, Rcpp::Named("iter") = iter,
		Rcpp::Named("info") = info);
}
