
#include "common.h"
#include "psoHeader.h"
#include "deHeader.h"
#include "psoFuncEval.h"

// RCPP FUNCTIONS
//[[Rcpp::export]]
Rcpp::List cppPSO(SEXP OBJ_FUNC, Rcpp::List PSO_INFO_LIST, const SEXP env, const bool IF_PARALLEL, const bool VERBOSE)
{
  //arma_rng::set_seed_random();
  /*int NCPU = omp_get_max_threads();
	if (NCPU < 3) NCPU = 2;
  omp_set_num_threads(NCPU - 1);*/

  // WRAP FUNCTIONS (INSPIRED BY R PACKAGE 'lbfgs')
  Rcpp::EvalBase *objfunc = NULL;
  Shield<SEXP> OBJ_FUNC_SEXP(Rcpp::as<SEXP>(OBJ_FUNC));
  if (TYPEOF(OBJ_FUNC_SEXP) == EXTPTRSXP) {
    objfunc = new Rcpp::EvalCompiled(OBJ_FUNC_SEXP, env);
  } else {
    objfunc = new Rcpp::EvalStandard(OBJ_FUNC_SEXP, env);
  }

  PSO_OPTIONS PSO_OPT; getAlgStruct(PSO_OPT, PSO_INFO_LIST);

  PSO_Result Result;
  if (VERBOSE) Rprintf("\nCalling Cpp PSO Kernel... ");
  PSO_MAIN(PSO_OPT, objfunc, IF_PARALLEL, VERBOSE, Result);
  if (VERBOSE) Rprintf("Done.\n");

  return List::create(Named("GBest") = wrap(Result.GBest),
                      Named("fGBest") = wrap(Result.fGBest),
                      Named("fGBestHist") = wrap(Result.fGBestHist),
                      Named("PBest") = wrap(Result.PBest),
                      Named("fPBest") = wrap(Result.fPBest));
}


// RCPP FUNCTIONS
//[[Rcpp::export]]
Rcpp::List cppDE(SEXP OBJ_FUNC, Rcpp::List DE_INFO_LIST, const SEXP env, const bool IF_PARALLEL, const bool VERBOSE)
{
  // arma_rng::set_seed_random();
  /*int NCPU = omp_get_max_threads();
   if (NCPU < 3) NCPU = 2;
   omp_set_num_threads(NCPU - 1);*/
  
  // WRAP FUNCTIONS (INSPIRED BY R PACKAGE 'lbfgs')
  Rcpp::EvalBase *objfunc = NULL;
  Shield<SEXP> OBJ_FUNC_SEXP(Rcpp::as<SEXP>(OBJ_FUNC));
  if (TYPEOF(OBJ_FUNC_SEXP) == EXTPTRSXP) {
    objfunc = new Rcpp::EvalCompiled(OBJ_FUNC_SEXP, env);
  } else {
    objfunc = new Rcpp::EvalStandard(OBJ_FUNC_SEXP, env);
  }
  
  DE_OPTIONS DE_OPT; getParamDE(DE_OPT, DE_INFO_LIST);

  DE_Result Result;
  if (VERBOSE) Rprintf("\nCalling Cpp DE Kernel... ");
  DE_MAIN(DE_OPT, objfunc, IF_PARALLEL, VERBOSE, Result);
  if (VERBOSE) Rprintf("Done.\n");
  
  return List::create(Named("GBest") = wrap(Result.GBest),
                      Named("fGBest") = wrap(Result.fGBest),
                      Named("fGBestHist") = wrap(Result.fGBestHist),
                      Named("PBest") = wrap(Result.PBest),
                      Named("fPBest") = wrap(Result.fPBest));
}


