
#include "psoHeader.h"

// RCPP FUNCTIONS
//[[Rcpp::export]]
Rcpp::List cppPSO(const int LOOPID, Rcpp::List PSO_INFO_LIST, Rcpp::List LBFGS_INFO_LIST,
                  Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST,
                  Rcpp::List EXTERNAL_LIST, const SEXP env, const bool IF_PARALLEL, const bool VERBOSE)
{
  //arma_rng::set_seed_random();
  /*int NCPU = omp_get_max_threads();
	if (NCPU < 3) NCPU = 2;
  omp_set_num_threads(NCPU - 1);*/
  OBJ_INFO OBJ; getInfoStruct(OBJ, OBJ_INFO_LIST);

  // WRAP FUNCTIONS (INSPIRED BY R PACKAGE 'lbfgs')
  Rcpp::EvalBase *dfnc = NULL;
  Shield<SEXP> DIST_FUNC_SEXP(Rcpp::as<SEXP>(OBJ_INFO_LIST["dist_func"]));
  if (TYPEOF(DIST_FUNC_SEXP) == EXTPTRSXP) {
    dfnc = new Rcpp::EvalCompiled(DIST_FUNC_SEXP, env);
  } else {
    dfnc = new Rcpp::EvalStandard(DIST_FUNC_SEXP, env);
  }

  int N_PAIR = OBJ.N_PAIR;
  arma::imat MODEL_PAIR = OBJ.MODEL_PAIR;

  model_diff_func *model_diff_ptr[N_PAIR];
  model_diff_func model_diff_collect[N_PAIR];

  for (int i = 0; i < N_PAIR; i++) {
    int tmID = MODEL_PAIR(i, 0); int rmID = MODEL_PAIR(i, 1);
    Rcpp::EvalBase *m1 = NULL;
    Shield<SEXP> tmp1(Rcpp::as<SEXP>(MODEL_INFO_LIST[tmID]));
    if (TYPEOF(tmp1) == EXTPTRSXP) {
      m1 = new Rcpp::EvalCompiled(tmp1, env);
    } else {
      m1 = new Rcpp::EvalStandard(tmp1, env);
    }
    Rcpp::EvalBase *m2 = NULL;
    Shield<SEXP> tmp2(Rcpp::as<SEXP>(MODEL_INFO_LIST[rmID]));
    if (TYPEOF(tmp2) == EXTPTRSXP) {
      m2 = new Rcpp::EvalCompiled(tmp2, env);
    } else {
      m2 = new Rcpp::EvalStandard(tmp2, env);
    }
    model_diff_collect[i] = model_diff_func(m1, m2, dfnc);
    model_diff_ptr[i] = &model_diff_collect[i];
  }

  PSO_OPTIONS PSO_OPT[N_PSO_OPTS]; getAlgStruct(PSO_OPT, PSO_INFO_LIST);
  LBFGS_PARAM LBFGS_OPTION; getNewtonStruct(LBFGS_OPTION, LBFGS_INFO_LIST);

  best_alpha_info external = {};
  if (OBJ.d_type == 1001) {
    external.DESIGN = as<arma::mat>(EXTERNAL_LIST["DESIGN"]);
    external.CRIT_VAL = as<double>(EXTERNAL_LIST["CRIT_VAL"]);
  }

  PSO_Result Result;
  if (VERBOSE) Rprintf("\nCalling Cpp PSO Kernel... ");
  PSO_MAIN(LOOPID, PSO_OPT, LBFGS_OPTION, OBJ, model_diff_ptr, &external, IF_PARALLEL, VERBOSE, Result);
  if (VERBOSE) Rprintf("Done.\n");

  return List::create(Named("GBest") = wrap(Result.GBest),
                      Named("fGBest") = wrap(Result.fGBest),
                      Named("fGBestHist") = wrap(Result.fGBestHist),
                      Named("PBest") = wrap(Result.PBest),
                      Named("fPBest") = wrap(Result.fPBest));
}


//[[Rcpp::export]]
Rcpp::List cppDesignCriterion(Rcpp::List PSO_INFO_LIST, Rcpp::List LBFGS_INFO_LIST,
                              Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST,
                              Rcpp::List EXTERNAL_LIST, SEXP env, arma::rowvec DESIGN)
{
  OBJ_INFO OBJ; getInfoStruct(OBJ, OBJ_INFO_LIST);

  // WRAP FUNCTIONS (INSPIRED BY R PACKAGE 'lbfgs')
  Rcpp::EvalBase *dfnc = NULL;
  Shield<SEXP> DIST_FUNC_SEXP(Rcpp::as<SEXP>(OBJ_INFO_LIST["dist_func"]));
  if (TYPEOF(DIST_FUNC_SEXP) == EXTPTRSXP) {
    dfnc = new Rcpp::EvalCompiled(DIST_FUNC_SEXP, env);
  } else {
    dfnc = new Rcpp::EvalStandard(DIST_FUNC_SEXP, env);
  }

  int N_PAIR = OBJ.N_PAIR;
  arma::imat MODEL_PAIR = OBJ.MODEL_PAIR;

  model_diff_func *model_diff_ptr[N_PAIR];
  model_diff_func model_diff_collect[N_PAIR];

  for (int i = 0; i < N_PAIR; i++) {
    int tmID = MODEL_PAIR(i, 0); int rmID = MODEL_PAIR(i, 1);
    Rcpp::EvalBase *m1 = NULL;
    Shield<SEXP> tmp1(Rcpp::as<SEXP>(MODEL_INFO_LIST[tmID]));
    if (TYPEOF(tmp1) == EXTPTRSXP) {
      m1 = new Rcpp::EvalCompiled(tmp1, env);
    } else {
      m1 = new Rcpp::EvalStandard(tmp1, env);
    }
    Rcpp::EvalBase *m2 = NULL;
    Shield<SEXP> tmp2(Rcpp::as<SEXP>(MODEL_INFO_LIST[rmID]));
    if (TYPEOF(tmp2) == EXTPTRSXP) {
      m2 = new Rcpp::EvalCompiled(tmp2, env);
    } else {
      m2 = new Rcpp::EvalStandard(tmp2, env);
    }
    model_diff_collect[i] = model_diff_func(m1, m2, dfnc);
    model_diff_ptr[i] = &model_diff_collect[i];
  }

  PSO_OPTIONS PSO_OPT[N_PSO_OPTS]; getAlgStruct(PSO_OPT, PSO_INFO_LIST);
  LBFGS_PARAM LBFGS_OPTION; getNewtonStruct(LBFGS_OPTION, LBFGS_INFO_LIST);

  arma::mat R_PARA;
  double val = DesignCriterion(0, PSO_OPT, LBFGS_OPTION, OBJ, model_diff_ptr, NULL, DESIGN, R_PARA);

  return List::create(Named("val") = wrap(val),
                      Named("theta2") = wrap(R_PARA));
}

//[[Rcpp::export]]
Rcpp::List cppEquivalence(Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST,
                          const double GBEST_VAL, const arma::mat PARA_SET, const arma::rowvec alpha, const SEXP env, const int nGrid)
{
  OBJ_INFO OBJ; getInfoStruct(OBJ, OBJ_INFO_LIST);

  // WRAP FUNCTIONS (INSPIRED BY R PACKAGE 'lbfgs')
  Rcpp::EvalBase *dfnc = NULL;
  Shield<SEXP> DIST_FUNC_SEXP(Rcpp::as<SEXP>(OBJ_INFO_LIST["dist_func"]));
  if (TYPEOF(DIST_FUNC_SEXP) == EXTPTRSXP) {
    dfnc = new Rcpp::EvalCompiled(DIST_FUNC_SEXP, env);
  } else {
    dfnc = new Rcpp::EvalStandard(DIST_FUNC_SEXP, env);
  }

  int N_PAIR = OBJ.N_PAIR;
  arma::imat MODEL_PAIR = OBJ.MODEL_PAIR;

  model_diff_func *model_diff_ptr[N_PAIR];
  model_diff_func model_diff_collect[N_PAIR];

  for (int i = 0; i < N_PAIR; i++) {
    int tmID = MODEL_PAIR(i, 0); int rmID = MODEL_PAIR(i, 1);
    Rcpp::EvalBase *m1 = NULL;
    Shield<SEXP> tmp1(Rcpp::as<SEXP>(MODEL_INFO_LIST[tmID]));
    if (TYPEOF(tmp1) == EXTPTRSXP) {
      m1 = new Rcpp::EvalCompiled(tmp1, env);
    } else {
      m1 = new Rcpp::EvalStandard(tmp1, env);
    }
    Rcpp::EvalBase *m2 = NULL;
    Shield<SEXP> tmp2(Rcpp::as<SEXP>(MODEL_INFO_LIST[rmID]));
    if (TYPEOF(tmp2) == EXTPTRSXP) {
      m2 = new Rcpp::EvalCompiled(tmp2, env);
    } else {
      m2 = new Rcpp::EvalStandard(tmp2, env);
    }
    model_diff_collect[i] = model_diff_func(m1, m2, dfnc);
    model_diff_ptr[i] = &model_diff_collect[i];
  }

  /* REVISE HERE */
  arma::vec xLine_1 = linspace<vec>(OBJ.dsLower(0), OBJ.dsUpper(0), nGrid);
  arma::vec xLine_2(1, fill::zeros);

  arma::mat DISPVALS;
  int dSupp = OBJ.dSupp;
  if (dSupp == 1) {
    arma::mat dsGrid(nGrid, 1); dsGrid.col(0) = xLine_1;
    arma::rowvec DIV = directionalDerivative(OBJ, dsGrid, PARA_SET, alpha, model_diff_ptr);
    DISPVALS.set_size(1, nGrid);
    DISPVALS.row(0) = DIV - GBEST_VAL;
  } else if (dSupp == 2) {
    xLine_2.reset(); xLine_2 = linspace<vec>(OBJ.dsLower(1), OBJ.dsUpper(1), nGrid);
    DISPVALS.set_size(nGrid, nGrid);
    for (int i = 0; i < nGrid; i++) {
      arma::mat dsGrid(nGrid, 2); dsGrid.col(0).fill(xLine_1(i)); dsGrid.col(1) = xLine_2;
      arma::rowvec DIV = directionalDerivative(OBJ, dsGrid, PARA_SET, alpha, model_diff_ptr);
      DISPVALS.row(i) = DIV - GBEST_VAL;
    }
  }
  double MAX_DISP = DISPVALS.max();
  /* END REVISION */

  /* OUTPUT */
  return List::create(Named("Grid_1") = wrap(xLine_1),
                      Named("Grid_2") = wrap(xLine_2),
                      Named("DirDeriv") = wrap(DISPVALS),
                      Named("MAX_DD") = wrap(MAX_DISP));
}

//[[Rcpp::export]]
Rcpp::List cppFedorovWynn(Rcpp::List FED_INFO_LIST, Rcpp::List LBFGS_INFO_LIST, Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST,
                          const SEXP env, const bool VERBOSE)
{
  OBJ_INFO OBJ; getInfoStruct(OBJ, OBJ_INFO_LIST);

  // WRAP FUNCTIONS (INSPIRED BY R PACKAGE 'lbfgs')
  Rcpp::EvalBase *dfnc = NULL;
  Shield<SEXP> DIST_FUNC_SEXP(Rcpp::as<SEXP>(OBJ_INFO_LIST["dist_func"]));
  if (TYPEOF(DIST_FUNC_SEXP) == EXTPTRSXP) {
    dfnc = new Rcpp::EvalCompiled(DIST_FUNC_SEXP, env);
  } else {
    dfnc = new Rcpp::EvalStandard(DIST_FUNC_SEXP, env);
  }

  int N_PAIR = OBJ.N_PAIR;
  arma::imat MODEL_PAIR = OBJ.MODEL_PAIR;

  model_diff_func *model_diff_ptr[N_PAIR];
  model_diff_func model_diff_collect[N_PAIR];

  for (int i = 0; i < N_PAIR; i++) {
    int tmID = MODEL_PAIR(i, 0); int rmID = MODEL_PAIR(i, 1);
    Rcpp::EvalBase *m1 = NULL;
    Shield<SEXP> tmp1(Rcpp::as<SEXP>(MODEL_INFO_LIST[tmID]));
    if (TYPEOF(tmp1) == EXTPTRSXP) {
      m1 = new Rcpp::EvalCompiled(tmp1, env);
    } else {
      m1 = new Rcpp::EvalStandard(tmp1, env);
    }
    Rcpp::EvalBase *m2 = NULL;
    Shield<SEXP> tmp2(Rcpp::as<SEXP>(MODEL_INFO_LIST[rmID]));
    if (TYPEOF(tmp2) == EXTPTRSXP) {
      m2 = new Rcpp::EvalCompiled(tmp2, env);
    } else {
      m2 = new Rcpp::EvalStandard(tmp2, env);
    }
    model_diff_collect[i] = model_diff_func(m1, m2, dfnc);
    model_diff_ptr[i] = &model_diff_collect[i];
  }

  FED_PARAM FED_OPTION; getFedorovStruct(FED_OPTION, FED_INFO_LIST);
  LBFGS_PARAM LBFGS_OPTION; getNewtonStruct(LBFGS_OPTION, LBFGS_INFO_LIST);

  FED_Result Result;
  if (VERBOSE) Rprintf("\nCalling Cpp Fedorov-Wynn Kernel... ");
  FEDOROVWYNN_MAIN(FED_OPTION, LBFGS_OPTION, OBJ, model_diff_ptr, VERBOSE, Result);
  if (VERBOSE) Rprintf("Done.\n");

  return List::create(Named("DESIGN") = wrap(Result.DESIGN),
                      Named("WT") = wrap(Result.WT),
                      Named("F_VAL") = wrap(Result.F_VAL),
                      Named("R_PARA") = wrap(Result.R_PARA),
                      Named("fvalHist") = wrap(Result.fvalHist));
}

//[[Rcpp::export]]
Rcpp::List cppUnifApprox(Rcpp::List REMES_INFO_LIST, Rcpp::List LBFGS_INFO_LIST, Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST,
                         const SEXP env, const bool VERBOSE)
{
  OBJ_INFO OBJ; getInfoStruct(OBJ, OBJ_INFO_LIST);

  // WRAP FUNCTIONS (INSPIRED BY R PACKAGE 'lbfgs')
  Rcpp::EvalBase *dfnc = NULL;
  Shield<SEXP> DIST_FUNC_SEXP(Rcpp::as<SEXP>(OBJ_INFO_LIST["dist_func"]));
  if (TYPEOF(DIST_FUNC_SEXP) == EXTPTRSXP) {
    dfnc = new Rcpp::EvalCompiled(DIST_FUNC_SEXP, env);
  } else {
    dfnc = new Rcpp::EvalStandard(DIST_FUNC_SEXP, env);
  }

  int N_PAIR = OBJ.N_PAIR;
  arma::imat MODEL_PAIR = OBJ.MODEL_PAIR;

  model_diff_func *model_diff_ptr[N_PAIR];
  model_diff_func model_diff_collect[N_PAIR];

  for (int i = 0; i < N_PAIR; i++) {
    int tmID = MODEL_PAIR(i, 0); int rmID = MODEL_PAIR(i, 1);
    Rcpp::EvalBase *m1 = NULL;
    Shield<SEXP> tmp1(Rcpp::as<SEXP>(MODEL_INFO_LIST[tmID]));
    if (TYPEOF(tmp1) == EXTPTRSXP) {
      m1 = new Rcpp::EvalCompiled(tmp1, env);
    } else {
      m1 = new Rcpp::EvalStandard(tmp1, env);
    }
    Rcpp::EvalBase *m2 = NULL;
    Shield<SEXP> tmp2(Rcpp::as<SEXP>(MODEL_INFO_LIST[rmID]));
    if (TYPEOF(tmp2) == EXTPTRSXP) {
      m2 = new Rcpp::EvalCompiled(tmp2, env);
    } else {
      m2 = new Rcpp::EvalStandard(tmp2, env);
    }
    model_diff_collect[i] = model_diff_func(m1, m2, dfnc);
    model_diff_ptr[i] = &model_diff_collect[i];
  }

  REMES_PARAM REMES_OPTION; getRemesStruct(REMES_OPTION, REMES_INFO_LIST);
  LBFGS_PARAM LBFGS_OPTION; getNewtonStruct(LBFGS_OPTION, LBFGS_INFO_LIST);

  REMES_Result Result;
  if (VERBOSE) Rprintf("\nCalling Cpp Remes Kernel... ");
  REMES_MAIN(REMES_OPTION, LBFGS_OPTION, OBJ, model_diff_ptr, VERBOSE, Result);
  if (VERBOSE) Rprintf("Done.\n");

  return List::create(Named("DESIGN") = wrap(Result.DESIGN),
                      Named("DD_DEV") = wrap(Result.DD_DEV),
                      Named("R_PARA") = wrap(Result.R_PARA),
                      Named("CPolyVal") = wrap(Result.CPolyVal));
}

