
// BODY
// PSO MAIN FUNCTIONS
void DE_MAIN(DE_OPTIONS DE_OPTS, Rcpp::EvalBase *objfunc,
             const bool IF_PARALLEL, const bool COUNTER_ON, DE_Result &DE_Result)
{
	/* -- BEGIN -- */
  // GET PSO PARAMETERS
	int nPop    = DE_OPTS.nPop;
	int dPop    = DE_OPTS.dPop;
	int maxIter    = DE_OPTS.maxIter; 
	int hasInitPop = DE_OPTS.hasInitPop;
	double freeRun   = DE_OPTS.freeRun; 
	double tol       = DE_OPTS.tol; 
  arma::rowvec varUpper  = DE_OPTS.varUpper;
  arma::rowvec varLower  = DE_OPTS.varLower;

	// DECLARE VARIABLES
  arma::mat pop(nPop, dPop), newpop(nPop, dPop), vMuta(nPop, dPop), PBest(nPop, dPop);
  arma::rowvec GBest(dPop);
  arma::vec fpop(nPop), fnewpop(nPop), fPBest(nPop);//
  double fGBest;
  arma::uword GBestIdx;
  arma::rowvec fGBestHist(maxIter + 1, fill::zeros);
	DE_DYN DE_DYN;
  
  /* -- START INITIALIZATION -- */
  if (COUNTER_ON) { Rprintf("DE Loop: Initializing .. "); }
  // INITIALIZE RANDOM SWARM
  pop = arma::randu(nPop, dPop) % repmat(varUpper - varLower, nPop, 1) + repmat(varLower, nPop, 1);
  if (hasInitPop > 0) {
    arma::mat initPop = DE_OPTS.initPop;
    for (arma::uword i = 0; i < initPop.n_rows; i++) {
      pop.row(i) = initPop.row(i);
    }
  }
  arma::rowvec fixedDims  = DE_OPTS.fixedDims;
  for (int j = 0; j < dPop; j ++) {
    if (!std::isnan(fixedDims(j))) {
      pop.col(j).fill(fixedDims(j));
    }
  }
  // INITIALIZE VELOCITY
  vMuta.fill(0);
  // INITIALIZE OBJECTIVE FUNCTION VALUES
  psoFuncEval(IF_PARALLEL, objfunc, pop, fpop); 
  // INITIALIZE LOCAL BEST
  fPBest = fpop;	PBest = pop;
  // INITIALIZE GLOBAL BEST
  fGBest = fPBest.min(GBestIdx); 
	GBest = PBest.row(GBestIdx);	
  // INITIALIZE PSO DYNAMIC PARAMETERS
  deUpdateDynPara(DE_OPTS, -1, DE_DYN, pop, PBest, GBest, fpop, fPBest, fGBest);
	// SAVE INITIAL GLOBAL BEST VALUE
	fGBestHist(0) = fGBest;
	  // SET ITERATION COUNTER
  int t; 
  if (COUNTER_ON) Rprintf("OK \n"); 
  /* -- FINISH INITIALIZATION -- */

  /* -- START PSO LOOP -- */
  for (t = 0; t < maxIter; t++) {
  	// PRINT OUT PROGRESS
    if (COUNTER_ON) {
      if (t == 0)  Rprintf("DE Loop: Updating ..    "); 
      Rprintf("\b\b\b%2.0f%%", (double)((t+1)*100/maxIter)); 
      if (t == (maxIter - 1)) Rprintf("\n"); 
    }
    // UPDATE MUTATION OPERATION
		deUpdatePop(DE_OPTS, DE_DYN, fpop, PBest, GBest, varUpper, varLower, pop, vMuta, newpop);
    // UPDATE POPULATION
    deCheckPop(DE_OPTS, DE_DYN, varUpper, varLower, newpop);	
    // UPDATE OBJECTIVE FUNCTION VALUES
    psoFuncEval(IF_PARALLEL, objfunc, newpop, fnewpop);
    // UPDATE POPULATION
    if (any(fnewpop < fpop)) {
      arma::uvec RowChange = find(fnewpop < fpop);
      fpop.elem(RowChange) = fnewpop.elem(RowChange);
      pop.rows(RowChange) = newpop.rows(RowChange);
    }
    // UPDATE THE LOCAL AND GLOBAL BEST
    if (any(fpop < fPBest)) {
      arma::uvec RowChange = find(fpop < fPBest);
      fPBest.elem(RowChange) = fpop.elem(RowChange);
      PBest.rows(RowChange) = pop.rows(RowChange);
    }
    if (min(fPBest) < fGBest) {
      fGBest = fPBest.min(GBestIdx); GBest = PBest.row(GBestIdx); //PSO_DYN.succ_GB = 1;
    }		
    // UPDATE DE DYNAMIC PARAMETERS
    deUpdateDynPara(DE_OPTS, t, DE_DYN, pop, PBest, GBest, fpop, fPBest, fGBest);
    // SAVE CURRENT GLOBAL BEST VALUE
    fGBestHist(t+1) = fGBest; 
    // CHECK STOPPING CRITERION
    if (t > (int)(freeRun*maxIter)) { 
      if (std::abs(fGBest - fGBestHist(t)) < tol) { 
        fGBestHist.subvec(t+1, maxIter).fill(fGBest); t = maxIter; 
        if (COUNTER_ON) Rprintf(" The updating procedure converges. \n");
      }
    }
  }
  /* -- FINISH PSO LOOP -- */
 
 	/* -- OUTPUT -- */
  DE_Result.GBest = GBest;
  DE_Result.fGBest = fGBest;
  DE_Result.fGBestHist = fGBestHist;
  DE_Result.PBest = PBest;
  DE_Result.fPBest = fPBest;
  /* -- END -- */
}
