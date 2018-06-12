// DECLARE FUNCTIONS
#include "DesignCriterion.h"

// BODY
void psoFuncEval(const bool IF_PARALLEL, const int LOOPID, PSO_OPTIONS PSO_OPTS[], const LBFGS_PARAM LBFGS_OPTION, 
								 const OBJ_INFO OBJ, const PSO_DYN PSO_DYN, model_diff_func *MODEL_COLLECTOR[], void *PSO_EXT, const mat swarm, vec &fSwarm)
{	
  int nSwarm = (int)swarm.n_rows;
 /* if (IF_PARALLEL) { 
		// PARALLEL LOOP (DOES NOT WORK FOR Rcpp::EvalBase*. MAY REPLACE OpenMP BY RcppParallel IN THE FUTURE.)
		int iParallel;
		#pragma omp parallel private(iParallel) 
		{
		#pragma omp for
			for (iParallel = 0; iParallel < nSwarm; iParallel++) {
				rowvec PARTICLE = conv_to<rowvec>::from(swarm.row(iParallel));
				// Optimal Design Criteria
				arma::mat R_PARA;
				fSwarm(iParallel) = DesignCriterion(LOOPID, PSO_OPTS, LBFGS_OPTION, OBJ, MODEL_COLLECTOR, PSO_EXT, PARTICLE, R_PARA);
			}
			#pragma omp barrier
		}
  } else {*/
		// NON-PARALLEL LOOP
		for (int iSwarm = 0; iSwarm < nSwarm; iSwarm++) {
			rowvec PARTICLE = arma::conv_to<rowvec>::from(swarm.row(iSwarm));
			// Optimal Design Criteria
			arma::mat R_PARA;
			fSwarm(iSwarm) = DesignCriterion(LOOPID, PSO_OPTS, LBFGS_OPTION, OBJ, MODEL_COLLECTOR, PSO_EXT, PARTICLE, R_PARA);
		}
  //}
}
