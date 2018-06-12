
// BODY
// FIND MINIMAL VALUE OF THE DIRECTIONAL DERIVATIVE FUNCTION
double max_CPoly_x(const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ, model_diff_func *MODEL_COLLECTOR[], const int PAIRID,
                   const arma::rowvec R_PARA, const arma::mat DESIGN, arma::rowvec &X_OUT)
{
  model_diff_func* func_input = MODEL_COLLECTOR[PAIRID];

  arma::imat MODEL_PAIR = OBJ.MODEL_PAIR;
  int tmID = MODEL_PAIR(PAIRID, 0);
  arma::rowvec T_PARA = OBJ.paras.submat(tmID, 0, tmID, OBJ.dParas(tmID) - 1);

  int dSupp = OBJ.dSupp;
  arma::rowvec X_UPPER = OBJ.dsUpper;
  arma::rowvec X_LOWER = OBJ.dsLower;
  arma::irowvec X_NBD(dSupp); X_NBD.fill(2);

  lbfgs_eval LBFGS_EVAL = {};
  LBFGS_EVAL.func_input = func_input;
  LBFGS_EVAL.UPPER    = X_UPPER;
  LBFGS_EVAL.LOWER    = X_LOWER;
  LBFGS_EVAL.NBD      = X_NBD;
  LBFGS_EVAL.T_PARA     = T_PARA;
  LBFGS_EVAL.R_PARA     = R_PARA;
  LBFGS_EVAL.FD_DELTA   = LBFGS_OPTION.FD_DELTA;

  lbfgs_parameter_t LBFGS_PAR;
  lbfgs_parameter_init(&LBFGS_PAR);
  LBFGS_PAR.m               = LBFGS_OPTION.LBFGS_LM;
  LBFGS_PAR.max_iterations  = LBFGS_OPTION.LBFGS_MAXIT;
  LBFGS_PAR.max_linesearch  = LBFGS_OPTION.LINESEARCH_MAXTRIAL;
  LBFGS_PAR.epsilon         = (lbfgsfloatval_t)LBFGS_OPTION.GRAD_EPS;
  LBFGS_PAR.delta           = (lbfgsfloatval_t)LBFGS_OPTION.FVAL_EPS;
  LBFGS_PAR.gtol            = (lbfgsfloatval_t)LBFGS_OPTION.LINESEARCH_WOLFE;
  LBFGS_PAR.ftol            = (lbfgsfloatval_t)LBFGS_OPTION.LINESEARCH_ARMIJO;
  LBFGS_PAR.min_step        = (lbfgsfloatval_t)LBFGS_OPTION.LINESEARCH_MIN;
  LBFGS_PAR.max_step        = (lbfgsfloatval_t)LBFGS_OPTION.LINESEARCH_MAX;

  lbfgsfloatval_t *X_VEC   = lbfgs_malloc(dSupp);
  lbfgsfloatval_t *X_VEC1  = lbfgs_malloc(dSupp);
  lbfgsfloatval_t fx, fx1; //int CONV;

  arma::rowvec X_INI = X_LOWER;
  for (int d = 0; d < dSupp; d++) { X_VEC[d] = (lbfgsfloatval_t)domainMapping(0, X_INI(d), X_NBD(d), X_UPPER(d), X_LOWER(d)); }
  lbfgs(dSupp, X_VEC, &fx, evaluate_remes_x, NULL, &LBFGS_EVAL, &LBFGS_PAR);

  //Rprintf("CONV: %d, f = %4.4f\n", CONV, fx);

  uword nSupp = DESIGN.n_rows;
  for (uword i = 0; i < (nSupp + 1); i++) {
    if (i < nSupp) { X_INI = DESIGN.row(i); } else { X_INI = X_UPPER; }
    for (int d = 0; d < dSupp; d++) { X_VEC1[d] = (lbfgsfloatval_t)domainMapping(0, X_INI(d), X_NBD(d), X_UPPER(d), X_LOWER(d)); }
    lbfgs(dSupp, X_VEC1, &fx1, evaluate_remes_x, NULL, &LBFGS_EVAL, &LBFGS_PAR);
    //Rprintf("CONV: %d, f = %4.4f\n", CONV, fx1);
    if (std::isfinite(fx1) & (!std::isnan(fx1))) {
      if ((!std::isfinite(fx)) | std::isnan(fx)) {
        fx = fx1; for (int d = 0; d < dSupp; d++) { X_VEC[d] = X_VEC1[d]; }
      } else {
        if (fx1 <= fx) {
          fx = fx1; for (int d = 0; d < dSupp; d++) { X_VEC[d] = X_VEC1[d]; }
        }
      }
    }
  }

  for (int d = 0; d < dSupp; d++) { X_OUT(d) = domainMapping(1, (double)X_VEC[d], X_NBD(d), X_UPPER(d), X_LOWER(d)); }

  lbfgs_free(X_VEC); lbfgs_free(X_VEC1);
  return (double)fx;
}

double min_CPoly_r(const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ, model_diff_func *MODEL_COLLECTOR[], const int PAIRID,
                   const arma::mat DESIGN, arma::rowvec &R_PARA_OUT)
{
  model_diff_func* func_input = MODEL_COLLECTOR[PAIRID];

  arma::imat MODEL_PAIR = OBJ.MODEL_PAIR;
  int tmID = MODEL_PAIR(PAIRID, 0);
  int rmID = MODEL_PAIR(PAIRID, 1);

  arma::rowvec T_PARA = OBJ.paras.submat(tmID, 0, tmID, OBJ.dParas(tmID) - 1);

  int dParas = OBJ.dParas(rmID);
  arma::rowvec R_UPPER = OBJ.parasUpper.submat(rmID, 0, rmID, dParas - 1);
  arma::rowvec R_LOWER = OBJ.parasLower.submat(rmID, 0, rmID, dParas - 1);
  arma::irowvec R_NBD(dParas);
  for (int d = 0; d < dParas; d++) { R_NBD(d) = (int)OBJ.parasBdd(rmID, d); }

  lbfgs_eval LBFGS_EVAL = {};
  LBFGS_EVAL.func_input = func_input;
  LBFGS_EVAL.UPPER    = OBJ.parasUpper.submat(rmID, 0, rmID, dParas - 1);
  LBFGS_EVAL.LOWER    = OBJ.parasLower.submat(rmID, 0, rmID, dParas - 1);
  LBFGS_EVAL.NBD      = R_NBD;
  LBFGS_EVAL.T_PARA     = T_PARA;
  LBFGS_EVAL.DESIGN     = DESIGN;
  LBFGS_EVAL.FD_DELTA   = LBFGS_OPTION.FD_DELTA;

  lbfgs_parameter_t LBFGS_PAR;
  lbfgs_parameter_init(&LBFGS_PAR);
  LBFGS_PAR.m               = LBFGS_OPTION.LBFGS_LM;
  LBFGS_PAR.max_iterations  = LBFGS_OPTION.LBFGS_MAXIT;
  LBFGS_PAR.max_linesearch  = LBFGS_OPTION.LINESEARCH_MAXTRIAL;
  LBFGS_PAR.epsilon         = (lbfgsfloatval_t)LBFGS_OPTION.GRAD_EPS;
  LBFGS_PAR.delta           = (lbfgsfloatval_t)LBFGS_OPTION.FVAL_EPS;
  LBFGS_PAR.gtol            = (lbfgsfloatval_t)LBFGS_OPTION.LINESEARCH_WOLFE;
  LBFGS_PAR.ftol            = (lbfgsfloatval_t)LBFGS_OPTION.LINESEARCH_ARMIJO;
  LBFGS_PAR.min_step        = (lbfgsfloatval_t)LBFGS_OPTION.LINESEARCH_MIN;
  LBFGS_PAR.max_step        = (lbfgsfloatval_t)LBFGS_OPTION.LINESEARCH_MAX;

  arma::rowvec R_PARA_INI = OBJ.parasInit.submat(rmID, 0, rmID, dParas - 1);

  lbfgsfloatval_t *R_PARA   = lbfgs_malloc(dParas);
  lbfgsfloatval_t *R_PARA1  = lbfgs_malloc(dParas);
  for (int d = 0; d < dParas; d++) { R_PARA[d] = (lbfgsfloatval_t)domainMapping(0, R_PARA_INI(d), R_NBD(d), R_UPPER(d), R_LOWER(d)); }

  lbfgsfloatval_t fx, fx1;

  int CONV;
  CONV = lbfgs(dParas, R_PARA, &fx, evaluate_remes_r, NULL, &LBFGS_EVAL, &LBFGS_PAR);
  //Rprintf("CONV: %d, f = %4.9f\n", CONV, fx);

  for (int d = 0; d < dParas; d++) { R_PARA1[d] = R_PARA[d]; }

  int LBFGS_RETRY = LBFGS_OPTION.LBFGS_RETRY; if (CONV) { LBFGS_RETRY++; }
  if ((!std::isfinite(fx)) | std::isnan(fx)) { LBFGS_RETRY++; }
  int count = 0;
  while (count < (2*LBFGS_RETRY)) {
    arma::rowvec SIGN_RAND = randu(1, dParas)*2.0 - 1.0;
    for (int d = 0; d < dParas; d++) { R_PARA1[d] = SIGN_RAND(d)*R_PARA1[d]; }
    CONV = lbfgs(dParas, R_PARA1, &fx1, evaluate_remes_r, NULL, &LBFGS_EVAL, &LBFGS_PAR);
    
    //Rprintf("CONV: %d, f = %4.9f\n", CONV, fx1);

    if (std::isfinite(fx1) & (!std::isnan(fx1))) {
      if ((!std::isfinite(fx)) | std::isnan(fx)) {
        fx = fx1;
        for (int d = 0; d < dParas; d++) { R_PARA[d] = R_PARA1[d]; }
      } else {
        if (fx1 < fx) {
          fx = fx1;
          for (int d = 0; d < dParas; d++) { R_PARA[d] = R_PARA1[d]; }
        }
      }
    }
    count++;
  }

  for (int d = 0; d < dParas; d++) { R_PARA_OUT(d) = domainMapping(1, (double)R_PARA[d], R_NBD(d), R_UPPER(d), R_LOWER(d)); }

  lbfgs_free(R_PARA); lbfgs_free(R_PARA1);
  return (double)fx;
}

double getCPolyVal(const OBJ_INFO OBJ, model_diff_func *MODEL_COLLECTOR[], const int PAIRID, const arma::mat DESIGN, const arma::rowvec R_PARA)
{
  model_diff_func* func_input = MODEL_COLLECTOR[PAIRID];

  arma::imat MODEL_PAIR = OBJ.MODEL_PAIR;
  int tmID = MODEL_PAIR(PAIRID, 0);
  arma::rowvec T_PARA = OBJ.paras.submat(tmID, 0, tmID, OBJ.dParas(tmID) - 1);

  return cpoly(R_PARA, T_PARA, DESIGN, func_input);
}

arma::rowvec getRivalDev(const OBJ_INFO OBJ, model_diff_func *MODEL_COLLECTOR[], const int PAIRID,
                         const arma::mat DESIGN, const arma::rowvec R_PARA)
{
  model_diff_func* func_input = MODEL_COLLECTOR[PAIRID];

  arma::imat MODEL_PAIR = OBJ.MODEL_PAIR;
  int rmID = MODEL_PAIR(PAIRID, 1);
  int tmID = MODEL_PAIR(PAIRID, 0);
  arma::rowvec T_PARA = OBJ.paras.submat(tmID, 0, tmID, OBJ.dParas(tmID) - 1);

  int dParas = OBJ.dParas(rmID);
  arma::rowvec R_UPPER = OBJ.parasUpper.submat(rmID, 0, rmID, dParas - 1);
  arma::rowvec R_LOWER = OBJ.parasLower.submat(rmID, 0, rmID, dParas - 1);
  arma::irowvec BDD(dParas, fill::zeros);

  double DELTA = 1e-3;
  for (int i = 0; i < dParas; i++) {
    if ((R_PARA(i) - 0.5*DELTA) < R_LOWER(i)) {
      BDD(i) = 1;
    } else if ((R_PARA(i) + 0.5*DELTA) > R_UPPER(i)) {
      BDD(i) = 2;
    }
  }

  double mul = cpoly(R_PARA, T_PARA, DESIGN, func_input);
  arma::rowvec rDev = rival_gr(R_PARA, DESIGN, func_input, BDD, DELTA);

  return mul*rDev;
}

// REMES MAIN FUNCTION
void REMES_MAIN(const REMES_PARAM REMES_OPTION, const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ, model_diff_func *MODEL_COLLECTOR[],
                const bool COUNTER_ON, REMES_Result &REMES_Result)
{
	/* -- BEGIN -- */
  // GET FEDOROV PARAMETERS
  int REMES_MAXIT = REMES_OPTION.REMES_MAXIT;
  // -- STOPPING CRITERION SETTING
  double freeRun = REMES_OPTION.freeRun;
  double REMES_EPS = REMES_OPTION.REMES_EPS;

  /* -- START INITIALIZATION -- */
  int nSupp = OBJ.nSupp;
  int dSupp = OBJ.dSupp;
  arma::mat DESIGN(nSupp, dSupp, fill::zeros);
  DESIGN = arma::randu(nSupp, dSupp) % arma::repmat(OBJ.dsUpper - OBJ.dsLower, nSupp, 1) + arma::repmat(OBJ.dsLower, nSupp, 1);
  DESIGN = arma::sort(DESIGN);
  DESIGN.row(0) = OBJ.dsLower;
  DESIGN.row(nSupp - 1) = OBJ.dsUpper;

  arma::rowvec R_PARA(OBJ.dParas(1));
  // INITIALIZE OBJECTIVE FUNCTION VALUES
  min_CPoly_r(LBFGS_OPTION, OBJ, MODEL_COLLECTOR, 0, DESIGN, R_PARA);
  //
	arma::mat ONE_PT(1, dSupp);
  arma::rowvec CPolyVal(nSupp);
	// SET ITERATION COUNTER
  int t;
  /* -- FINISH INITIALIZATION -- */

  /* -- START REMES LOOP -- */
  for (t = 0; t < REMES_MAXIT; t++) {
  	// PRINT OUT PROGRESS
    if (COUNTER_ON) {
      if (t == 0)  Rprintf("Remes: Updating ..    ");
      Rprintf("\b\b\b%2.0f%%", (double)((t+1)*100/REMES_MAXIT));
      if (t == (REMES_MAXIT - 1)) Rprintf("\n");
    }

    arma::rowvec X_0(dSupp);
    max_CPoly_x(LBFGS_OPTION, OBJ, MODEL_COLLECTOR, 0, R_PARA, DESIGN, X_0);
    //
    // Rprintf("BEFORE:\n"); matrixPrintf(DESIGN); rvecPrintf(X_0);
    double s0, s1;
    ONE_PT.row(0) = X_0;
    s0 = getCPolyVal(OBJ, MODEL_COLLECTOR, 0, ONE_PT, R_PARA);
    if (X_0[0] < DESIGN(0, 0)) {
      // x0 is at the most left position
      ONE_PT.row(0) = DESIGN.row(0);
      s1 = getCPolyVal(OBJ, MODEL_COLLECTOR, 0, ONE_PT, R_PARA);
      if (s0*s1 > 0) { DESIGN.row(0) = X_0; } else { DESIGN.row(nSupp - 1) = X_0; }
    } else if (X_0[0] >= DESIGN(nSupp - 1, 0)) {
      // x0 is at the most right position
      ONE_PT.row(0) = DESIGN.row(nSupp - 1);
      s1 = getCPolyVal(OBJ, MODEL_COLLECTOR, 0, ONE_PT, R_PARA);
      if (s0*s1 > 0) { DESIGN.row(nSupp - 1) = X_0; } else { DESIGN.row(0) = X_0; }
    } else {
      for (int i = 0; i < (nSupp - 1); i++) {
        if ((X_0[0] >= DESIGN(i, 0)) & (X_0[0] < DESIGN(i+1, 0))) {
          ONE_PT.row(0) = DESIGN.row(i);
          s1 = getCPolyVal(OBJ, MODEL_COLLECTOR, 0, ONE_PT, R_PARA);
          if (s0*s1 > 0) { DESIGN.row(i) = X_0; } else { DESIGN.row(i+1) = X_0; }
          i = nSupp;
        }
      }
    }
    DESIGN = arma::sort(DESIGN); // Rprintf("AFTER:\n");matrixPrintf(DESIGN);

    min_CPoly_r(LBFGS_OPTION, OBJ, MODEL_COLLECTOR, 0, DESIGN, R_PARA);

    // CHECK STOPPING CRITERION
    if (t > (int)(freeRun*REMES_MAXIT)) {
      for (int i = 0; i < nSupp; i++) {
        ONE_PT.row(0) = DESIGN.row(i);
        CPolyVal(i) = getCPolyVal(OBJ, MODEL_COLLECTOR, 0, ONE_PT, R_PARA);
      }
      CPolyVal = arma::abs(CPolyVal);
      double REMES_CHECK = 1.0 - CPolyVal.min()/CPolyVal.max();
      if (REMES_CHECK < REMES_EPS) {
        t = REMES_MAXIT;
        if (COUNTER_ON) Rprintf(" The updating procedure converges. \n");
      }
    }
  }
  /* -- FINISH FEDOROV LOOP -- */
  // Calculate Design Weights by (2.8) in Dette and Titoff (2009)
  arma::mat DD_DEV(R_PARA.n_elem, nSupp, fill::zeros);
  for (int i = 0; i < nSupp; i++) {
    ONE_PT.row(0) = DESIGN.row(i);
    CPolyVal(i) = getCPolyVal(OBJ, MODEL_COLLECTOR, 0, ONE_PT, R_PARA);
    DD_DEV.submat(0, i, R_PARA.n_elem - 1, i) = getRivalDev(OBJ, MODEL_COLLECTOR, 0, ONE_PT, R_PARA).t();
  }
 	/* -- OUTPUT -- */
  REMES_Result.DESIGN = DESIGN;
  REMES_Result.DD_DEV = DD_DEV;
  REMES_Result.R_PARA = R_PARA;
  REMES_Result.CPolyVal = CPolyVal;
  /* -- END -- */
}

