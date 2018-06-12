
// BODY
// FIND MINIMAL VALUE OF THE DIRECTIONAL DERIVATIVE FUNCTION
double minDirDev(const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ, model_diff_func *MODEL_COLLECTOR[], const int PAIRID,
                 const arma::rowvec R_PARA, arma::rowvec &X_OUT)
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

  arma::rowvec X_INI = randu(1, dSupp) % (X_UPPER - X_LOWER) + X_LOWER;

  lbfgsfloatval_t *X_VEC   = lbfgs_malloc(dSupp);
  lbfgsfloatval_t *X_VEC1  = lbfgs_malloc(dSupp);
  for (int d = 0; d < dSupp; d++) { X_VEC[d] = (lbfgsfloatval_t)domainMapping(0, X_INI(d), X_NBD(d), X_UPPER(d), X_LOWER(d)); }

  lbfgsfloatval_t fx, fx1;

  int CONV;
  CONV = lbfgs(dSupp, X_VEC, &fx, evaluate_dirdev, NULL, &LBFGS_EVAL, &LBFGS_PAR);
  //Rprintf("CONV: %d\n", CONV);
  int LBFGS_RETRY = LBFGS_OPTION.LBFGS_RETRY; if (CONV) { LBFGS_RETRY++; }
  if ((!std::isfinite(fx)) | std::isnan(fx)) { LBFGS_RETRY++; }
  int count = 0;
  while ((count < LBFGS_RETRY)) {
    X_INI = randu(1, dSupp) % (X_UPPER - X_LOWER) + X_LOWER;
    for (int d = 0; d < dSupp; d++) { X_VEC1[d] = (lbfgsfloatval_t)domainMapping(0, X_INI(d), X_NBD(d), X_UPPER(d), X_LOWER(d)); }
    CONV = lbfgs(dSupp, X_VEC1, &fx1, evaluate_dirdev, NULL, &LBFGS_EVAL, &LBFGS_PAR);
    //Rprintf("CONV: %d\n", CONV);
    if (std::isfinite(fx1) & (!std::isnan(fx1))) {
      if ((!std::isfinite(fx)) | std::isnan(fx)) {
        fx = fx1;
        for (int d = 0; d < dSupp; d++) { X_VEC[d] = X_VEC1[d]; }
      } else {
        if (fx1 < fx) {
          fx = fx1;
          for (int d = 0; d < dSupp; d++) { X_VEC[d] = X_VEC1[d]; }
        }
      }
    }
    count++;
  }
  for (int d = 0; d < dSupp; d++) { X_OUT(d) = domainMapping(1, (double)X_VEC[d], X_NBD(d), X_UPPER(d), X_LOWER(d)); }

  lbfgs_free(X_VEC); lbfgs_free(X_VEC1);
  return (double)fx;
}

// Trimming
void trimDesign(arma::mat &DESIGN, arma::rowvec &WT, const double tol)
{
  uword n = DESIGN.n_rows;
  uword d = DESIGN.n_cols;
  arma::ivec GROUP_ID(n); GROUP_ID.fill(-1);
  arma::rowvec DIFF_TWO_PTS(d);
  
  int gid = 0;
  for (uword i = 0; i < n; i++) {
    if (GROUP_ID(i) < 0) {
      GROUP_ID(i) = gid;
      for (uword j = 0; j < n; j++) {
        if (GROUP_ID(j) < 0) {
          DIFF_TWO_PTS = arma::abs(DESIGN.row(i) - DESIGN.row(j));
          if (arma::max(DIFF_TWO_PTS) < tol) { GROUP_ID(j) = gid; }
        }
      }
      gid++;
    }
  }

  arma::mat NEW_DESIGN(gid, d, fill::zeros);
  arma::rowvec NEW_WT(gid, fill::zeros);
  arma::uvec g_i;
  for (int i = 0; i < gid; i++) {
    g_i.reset(); g_i = arma::find(GROUP_ID == i);
    for (uword j = 0; j < g_i.n_elem; j++) {
      NEW_DESIGN.row(i) = NEW_DESIGN.row(i) + DESIGN.row(g_i(j));
      NEW_WT(i) = NEW_WT(i) + WT(g_i(j));
    }
    double ng_double = (double)g_i.n_elem;
    NEW_DESIGN.row(i) = NEW_DESIGN.row(i)/ng_double;
  }

  DESIGN.set_size(arma::size(NEW_DESIGN)); DESIGN = NEW_DESIGN;
  WT.set_size(arma::size(NEW_WT)); WT = NEW_WT;
}

// FIND MINIMAL VALUE OF THE DIRECTIONAL DERIVATIVE FUNCTION
double optAlpha(const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ, model_diff_func *MODEL_COLLECTOR[],
                const arma::mat DESIGN, const arma::rowvec WT, const arma::vec ALPHA_CAND)
{
  uword n = ALPHA_CAND.n_elem;

  arma::rowvec ALPHA_VEC(arma::size(WT)); 
  arma::rowvec NEW_WT(arma::size(WT));
  arma::vec VAL_CAND(n);

  for (uword i = 0; i < n; i++) {
    ALPHA_VEC.zeros(); NEW_WT.zeros();
    ALPHA_VEC.subvec(0, DESIGN.n_rows - 2).fill(1.0 - ALPHA_CAND[i]);  
    ALPHA_VEC(DESIGN.n_rows - 1) = ALPHA_CAND[i];
    NEW_WT = WT % ALPHA_VEC;
    arma::rowvec R_PARA(OBJ.dParas(1));
    VAL_CAND(i) = minDistCalc(LBFGS_OPTION, OBJ, MODEL_COLLECTOR, 0, DESIGN, NEW_WT, R_PARA);
  }
  uword ai = index_max(VAL_CAND);
  return ALPHA_CAND(ai);
}


// FEDOROV-WYNN MAIN FUNCTION
void FEDOROVWYNN_MAIN(const FED_PARAM FED_OPTION, const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ, model_diff_func *MODEL_COLLECTOR[],
                      const bool COUNTER_ON, FED_Result &FED_Result)
{
	/* -- BEGIN -- */
  // GET FEDOROV PARAMETERS
  int FED_MAXIT = FED_OPTION.FED_MAXIT;
  int FED_TRIM = FED_OPTION.FED_TRIM;
  double FED_TRIM_EPS = FED_OPTION.FED_TRIM_EPS;
  // -- STOPPING CRITERION SETTING
  double freeRun = FED_OPTION.freeRun;
  double FED_EPS = FED_OPTION.FED_EPS;
  // -- UPDATING PARAMETER
  int FED_ALPHA_GRID = FED_OPTION.FED_ALPHA_GRID;
  arma::vec ALPHA_CAND = linspace(0.0, 1.0, FED_ALPHA_GRID);

  /* -- START INITIALIZATION -- */
  int nSupp = OBJ.nSupp;
  int dSupp = OBJ.dSupp;
  arma::mat DESIGN(nSupp, dSupp, fill::zeros);
  DESIGN = randu(nSupp, dSupp) % repmat(OBJ.dsUpper - OBJ.dsLower, nSupp, 1) + repmat(OBJ.dsLower, nSupp, 1);
  arma::rowvec WT(nSupp); WT.fill((1.0/(double)nSupp));
  
  arma::mat NEW_DESIGN;
  arma::rowvec NEW_WT, NEW_ALPHA;

  arma::rowvec R_PARA(OBJ.dParas(1));
  // INITIALIZE OBJECTIVE FUNCTION VALUES
  double fval = minDistCalc(LBFGS_OPTION, OBJ, MODEL_COLLECTOR, 0, DESIGN, WT, R_PARA);
	// SAVE INITIAL GLOBAL BEST VALUE
  arma::rowvec fvalHist(FED_MAXIT + 1);
	fvalHist(0) = fval;
  //
	// SET ITERATION COUNTER
  int t; int trim = 0; arma::uvec WT_MATTERS;
  /* -- FINISH INITIALIZATION -- */

  /* -- START FEDOROV LOOP -- */
  for (t = 0; t < FED_MAXIT; t++) {
  	// PRINT OUT PROGRESS
    if (COUNTER_ON) {
      if (t == 0)  Rprintf("Fedorov-Wynn: Updating ..    "); 
      Rprintf("\b\b\b%2.0f%%", (double)((t+1)*100/FED_MAXIT)); 
      if (t == (FED_MAXIT - 1)) Rprintf("\n"); 
    }
    
    trim++;

    arma::rowvec NEW_X(dSupp);
    minDirDev(LBFGS_OPTION, OBJ, MODEL_COLLECTOR, 0, R_PARA, NEW_X);
    //
    
    NEW_DESIGN.set_size(DESIGN.n_rows + 1, dSupp);
    NEW_DESIGN.rows(0, DESIGN.n_rows - 1) = DESIGN;
    NEW_DESIGN.row(DESIGN.n_rows) = NEW_X;

    NEW_WT.set_size(NEW_DESIGN.n_rows);  
    NEW_WT.subvec(0, NEW_DESIGN.n_rows - 2) = WT;
    NEW_WT(NEW_DESIGN.n_rows - 1) = 1.0;
    
    double ALPHA_VAL = optAlpha(LBFGS_OPTION, OBJ, MODEL_COLLECTOR, NEW_DESIGN, NEW_WT, ALPHA_CAND);

    NEW_ALPHA.set_size(NEW_DESIGN.n_rows);  
    NEW_ALPHA.subvec(0, NEW_DESIGN.n_rows - 2).fill(1.0 - ALPHA_VAL);
    NEW_ALPHA(NEW_DESIGN.n_rows - 1) = ALPHA_VAL;
    NEW_WT = NEW_WT % NEW_ALPHA;

    DESIGN.set_size(arma::size(NEW_DESIGN)); DESIGN = NEW_DESIGN;
    WT.set_size(arma::size(NEW_WT)); WT = NEW_WT;

    if (trim == FED_TRIM) {
      trimDesign(DESIGN, WT, FED_TRIM_EPS); trim = 0;
      if (arma::any(WT <= 1e-12)) {
        WT_MATTERS.reset(); WT_MATTERS = arma::find(WT > 1e-12);

        NEW_DESIGN.set_size(WT_MATTERS.n_elem, dSupp);
        NEW_DESIGN = DESIGN.rows(WT_MATTERS);

        NEW_WT.set_size(NEW_DESIGN.n_rows);         
        for (uword k = 0; k < WT_MATTERS.n_elem; k++) { NEW_WT(k) = WT(WT_MATTERS(k)); }
        
        DESIGN.set_size(arma::size(NEW_DESIGN)); DESIGN = NEW_DESIGN;
        WT.set_size(arma::size(NEW_WT)); WT = NEW_WT;  
      }
    }
      
    fval = minDistCalc(LBFGS_OPTION, OBJ, MODEL_COLLECTOR, 0, DESIGN, WT, R_PARA);
    // SAVE CURRENT GLOBAL BEST VALUE
    fvalHist(t+1) = fval; 
    // CHECK STOPPING CRITERION
    if (t > (int)(freeRun*FED_MAXIT)) { 
      if (std::abs(fval - fvalHist(t)) < FED_EPS) { 
        fvalHist.subvec(t+1, FED_MAXIT).fill(fval); t = FED_MAXIT; 
        if (COUNTER_ON) Rprintf(" The updating procedure converges. \n");
      }
    }
  }

  trimDesign(DESIGN, WT, FED_TRIM_EPS);
  /* -- FINISH FEDOROV LOOP -- */

 	/* -- OUTPUT -- */
  FED_Result.DESIGN = DESIGN;
  FED_Result.WT = WT;
  FED_Result.F_VAL = fval;
  FED_Result.R_PARA = R_PARA;
  FED_Result.fvalHist = fvalHist;
  /* -- END -- */
}

