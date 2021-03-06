// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cppPSO
Rcpp::List cppPSO(SEXP OBJ_FUNC, Rcpp::List PSO_INFO_LIST, const SEXP env, const bool IF_PARALLEL, const bool VERBOSE);
RcppExport SEXP _globpso_cppPSO(SEXP OBJ_FUNCSEXP, SEXP PSO_INFO_LISTSEXP, SEXP envSEXP, SEXP IF_PARALLELSEXP, SEXP VERBOSESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type OBJ_FUNC(OBJ_FUNCSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type PSO_INFO_LIST(PSO_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< const SEXP >::type env(envSEXP);
    Rcpp::traits::input_parameter< const bool >::type IF_PARALLEL(IF_PARALLELSEXP);
    Rcpp::traits::input_parameter< const bool >::type VERBOSE(VERBOSESEXP);
    rcpp_result_gen = Rcpp::wrap(cppPSO(OBJ_FUNC, PSO_INFO_LIST, env, IF_PARALLEL, VERBOSE));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_globpso_cppPSO", (DL_FUNC) &_globpso_cppPSO, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_globpso(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
