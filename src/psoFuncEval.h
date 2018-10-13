// DECLARE FUNCTIONS

// BODY
void psoFuncEval(const bool IF_PARALLEL, Rcpp::EvalBase *objfunc, const mat swarm, vec &fSwarm)
{	
  int nSwarm = (int)swarm.n_rows;
  double feval = 0;
 /* if (IF_PARALLEL) { 
		// PARALLEL LOOP (DOES NOT WORK FOR Rcpp::EvalBase*. MAY REPLACE OpenMP BY RcppParallel IN THE FUTURE.)
		int iParallel;
		#pragma omp parallel private(iParallel) 
		{
		#pragma omp for
			for (iParallel = 0; iParallel < nSwarm; iParallel++) {
				rowvec PARTICLE = arma::conv_to<rowvec>::from(swarm.row(iParallel));
				Shield<SEXP> PARTICLE_SEXP(Rcpp::wrap(PARTICLE));
				Rcpp::NumericVector PARTICLE_Rform = Rcpp::as<Rcpp::NumericVector>(PARTICLE_SEXP);
				feval = (double) objfunc->eval(PARTICLE_Rform);
				fSwarm(iSwarm) = feval;
			}
			#pragma omp barrier
		}
  } else {*/
		// NON-PARALLEL LOOP
		for (int iSwarm = 0; iSwarm < nSwarm; iSwarm++) {
			rowvec PARTICLE = arma::conv_to<rowvec>::from(swarm.row(iSwarm));
			Shield<SEXP> PARTICLE_SEXP(Rcpp::wrap(PARTICLE));
			//Rcpp::NumericVector PARTICLE_Rform = Rcpp::as<Rcpp::NumericVector>(PARTICLE_SEXP);
			feval = (double) objfunc->eval(PARTICLE_SEXP);
			fSwarm(iSwarm) = feval;
		}
  //}
}
