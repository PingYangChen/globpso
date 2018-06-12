
// INCLUDE HEADER FILES
#include <malloc.h>
#include "lbfgsKernel.h"
#include "lbfgsEval.h"

// DECLARE FUNCTIONS
double criterionList(const int LOOPID, PSO_OPTIONS PSO_OPTS[], const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ,
                     model_diff_func *MODEL_COLLECTOR[],
                     const arma::mat DESIGN, const arma::rowvec WT, arma::mat &R_PARA);
double minDistCalc(const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ, model_diff_func *MODEL_COLLECTOR[], const int PAIRID,
                   const arma::mat DESIGN, const arma::rowvec WT, arma::rowvec &R_PARA_OUT);
arma::rowvec directionalDerivative(const OBJ_INFO OBJ, const arma::mat dsGrid, const arma::mat PARA_SET, const arma::rowvec alpha,
                                   model_diff_func *MODEL_COLLECTOR[]);
arma::rowvec distCalc(const OBJ_INFO OBJ, const arma::mat x, const arma::mat PARA_SET,
                      model_diff_func *MODEL_COLLECTOR[], const int PAIRID);

// BODY
double DesignCriterion(const int LOOPID, PSO_OPTIONS PSO_OPTS[], const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ,
                       model_diff_func *MODEL_COLLECTOR[], void *PSO_EXT, const rowvec x, arma::mat &R_PARA)
{
	int nSupp = OBJ.nSupp;
	int dSupp = OBJ.dSupp;
	int d_type = OBJ.d_type;
	double minWt = OBJ.minWt;
	//Rprintf("get design\n");
  arma::rowvec swarm = x;

  double val = 1e10;
  if (LOOPID == 0) {
  	switch (d_type) {
  		// Exact Design Module
  		case 0:
  		{
        arma::mat DESIGN(nSupp, dSupp, fill::zeros);
        for (int i = 0; i < nSupp; i++) { DESIGN.row(i) = swarm.subvec(i*dSupp, (i+1)*dSupp - 1); }
        arma::rowvec WT(nSupp, fill::zeros);
  			WT.fill(1.0/(double)nSupp);
        val = criterionList(LOOPID, PSO_OPTS, LBFGS_OPTION, OBJ, MODEL_COLLECTOR, DESIGN, WT, R_PARA);
        val *= -1.0;
  			break;
  		}
  		// Approximation Design Module
  		case 1:
  		{
        arma::mat DESIGN(nSupp, dSupp, fill::zeros);
        for (int i = 0; i < nSupp; i++) { DESIGN.row(i) = swarm.subvec(i*dSupp, (i+1)*dSupp - 1); }
        arma::rowvec WT(nSupp);
  			// Get Design Weight
  			arma::rowvec ang = swarm.subvec(nSupp*dSupp, nSupp*dSupp + nSupp - 2);
  			arma::rowvec wcumsin(nSupp, fill::zeros), wcos(nSupp, fill::zeros);
  			arma::rowvec wsin = arma::sin(ang);
  			wcumsin(0) = 1.0;
  			for (int i = 1; i < nSupp; i++) { wcumsin(i) = wcumsin(i-1)*wsin(i-1); }
  			wcos(nSupp - 1) = 1.0; wcos.subvec(0, nSupp - 2) = arma::cos(ang);
  			WT = wcumsin % wcos;
  			WT = WT % WT;
  			if (arma::min(WT) > minWt) {
          val = criterionList(LOOPID, PSO_OPTS, LBFGS_OPTION, OBJ, MODEL_COLLECTOR, DESIGN, WT, R_PARA);
          val *= -1.0;
  			}
  			break;
  		}
      // Find Optimal Weight for Equivalence Theorem on Max-min Optimal Design
      case 1001:
      {
        best_alpha_info EXT = *(best_alpha_info*)(PSO_EXT);
        arma::mat DESIGN = EXT.DESIGN;
        double CRIT_VAL = EXT.CRIT_VAL;

        int n_model = (int)swarm.n_elem;
        arma::rowvec alpha(n_model + 1, fill::zeros);
        arma::rowvec ang = swarm;
        arma::rowvec wcumsin(n_model + 1, fill::zeros), wcos(n_model + 1, fill::zeros);
        arma::rowvec wsin = arma::sin(ang);
        wcumsin(0) = 1.0;
        for (int i = 1; i < (n_model + 1); i++) { wcumsin(i) = wcumsin(i-1)*wsin(i-1); }
        wcos(n_model) = 1.0; wcos.subvec(0, n_model - 1) = arma::cos(ang);
        alpha = wcumsin % wcos;
        alpha = alpha % alpha;
        arma::rowvec DIV = directionalDerivative(OBJ, DESIGN, OBJ.paras, alpha, MODEL_COLLECTOR);
        DIV -= CRIT_VAL;
        val = arma::accu(DIV % DIV);
        break;
      }
  	}
  } else {
    inner_pso_info EXT = *(inner_pso_info*)(PSO_EXT);
    int PAIRID = (int) EXT.PAIRID;
    arma::mat DESIGN = (arma::mat) EXT.DESIGN;
    arma::rowvec WT = (arma::rowvec) EXT.WT;

    arma::imat MODEL_PAIR = (arma::imat) OBJ.MODEL_PAIR;
    int tmID = MODEL_PAIR(PAIRID, 0);
    arma::rowvec T_PARA = OBJ.paras.submat(tmID, 0, tmID, OBJ.dParas(tmID) - 1);

    model_diff_func* func_input = MODEL_COLLECTOR[PAIRID];
    Rcpp::EvalBase *m1_func = (Rcpp::EvalBase *) func_input->M1_FUNC;
    Rcpp::EvalBase *m2_func = (Rcpp::EvalBase *) func_input->M2_FUNC;
    Rcpp::EvalBase *distFunc = (Rcpp::EvalBase *) func_input->DISTFUNC;

    Shield<SEXP> DESIGN_SEXP(Rcpp::wrap(DESIGN));
    Shield<SEXP> T_PARA_SEXP(Rcpp::wrap(T_PARA));
    Shield<SEXP> R_PARA_SEXP(Rcpp::wrap(swarm));

    Rcpp::NumericMatrix DESIGN_Rform = Rcpp::as<Rcpp::NumericMatrix>(DESIGN_SEXP);
    Rcpp::NumericVector T_PARA_Rform = Rcpp::as<Rcpp::NumericVector>(T_PARA_SEXP);
    Rcpp::NumericVector R_PARA_Rform = Rcpp::as<Rcpp::NumericVector>(R_PARA_SEXP);

    Rcpp::NumericVector eta_T_Rform((int)WT.n_elem), eta_R_Rform((int)WT.n_elem), DIV_Rform((int)WT.n_elem);

    eta_T_Rform = (Rcpp::NumericVector) m1_func->eval(DESIGN_Rform, T_PARA_Rform);
    if (Rcpp::all(Rcpp::is_finite(eta_T_Rform))) {
      eta_R_Rform = (Rcpp::NumericVector) m2_func->eval(DESIGN_Rform, R_PARA_Rform);
      if (Rcpp::all(Rcpp::is_finite(eta_R_Rform))) {
        DIV_Rform = (Rcpp::NumericVector) distFunc->eval(eta_T_Rform, eta_R_Rform);
        if (Rcpp::all(Rcpp::is_finite(DIV_Rform))) {
          arma::rowvec DIV(DIV_Rform.begin(), DIV_Rform.size());
          val = 0;
          for (uword i = 0; i < WT.n_elem; i++) { val += WT(i)*DIV(i); }
          //val = arma::accu(WT % DIV);
        }
      }
    }
  }
	return val;
}

//
double criterionList(const int LOOPID, PSO_OPTIONS PSO_OPTS[], const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ,
                     model_diff_func *MODEL_COLLECTOR[],
										 const arma::mat DESIGN, const arma::rowvec WT, arma::mat &R_PARA)
{
	int crit_type = OBJ.crit_type;
	int N_PAIR = OBJ.N_PAIR;
  int LBFGS = LBFGS_OPTION.IF_INNER_LBFGS;

  R_PARA.reset(); R_PARA.set_size(OBJ.dParas.n_elem, OBJ.dParas.max()); R_PARA.zeros();
	double val = 1e10;
	switch (crit_type) {
		case 0: // Fixed True
		{
      int rmID = OBJ.MODEL_PAIR(0, 1);
      R_PARA.submat(0, 0, 0, OBJ.dParas(0) - 1) = OBJ.paras.submat(0, 0, 0, OBJ.dParas(0) - 1);
			arma::rowvec R_PARA_tmp(OBJ.dParas(rmID));

      if (LBFGS == 0) {
        inner_pso_info PSO_EXT = {};
        PSO_EXT.PAIRID = 0;
        PSO_EXT.DESIGN = DESIGN;
        PSO_EXT.WT = WT;

        PSO_OPTS[LOOPID + 1].dSwarm = OBJ.dParas(rmID);
        PSO_OPTS[LOOPID + 1].varUpper.set_size(OBJ.dParas(rmID));
        PSO_OPTS[LOOPID + 1].varUpper = OBJ.parasUpper.submat(rmID, 0, rmID, OBJ.dParas(rmID) - 1);
        PSO_OPTS[LOOPID + 1].varLower.set_size(OBJ.dParas(rmID));
        PSO_OPTS[LOOPID + 1].varLower = OBJ.parasLower.submat(rmID, 0, rmID, OBJ.dParas(rmID) - 1);

        PSO_Result InnerResult;
        PSO_MAIN(LOOPID + 1, PSO_OPTS, LBFGS_OPTION, OBJ, MODEL_COLLECTOR, &PSO_EXT, FALSE, FALSE, InnerResult);
        R_PARA_tmp = InnerResult.GBest;
        val = InnerResult.fGBest;
      } else {
        val = minDistCalc(LBFGS_OPTION, OBJ, MODEL_COLLECTOR, 0, DESIGN, WT, R_PARA_tmp);
      }

      R_PARA.submat(1, 0, 1, OBJ.dParas(rmID) - 1) = R_PARA_tmp;
			break;
		}
		case 1: // Max-min, Fixed True
		{
      R_PARA.submat(0, 0, 0, OBJ.dParas(0) - 1) = OBJ.paras.submat(0, 0, 0, OBJ.dParas(0) - 1);
			arma::rowvec std_vals = OBJ.std_vals;
			arma::rowvec eff_vals(N_PAIR, fill::zeros);
			for (int i = 0; i < N_PAIR; i++) {
				arma::rowvec R_PARA_tmp(OBJ.dParas(i+1));

        if (LBFGS == 0) {
          inner_pso_info PSO_EXT = {};
          PSO_EXT.PAIRID = i;
          PSO_EXT.DESIGN = DESIGN;
          PSO_EXT.WT = WT;

          int rmID = OBJ.MODEL_PAIR(i, 1);

          PSO_OPTS[LOOPID + 1].dSwarm = OBJ.dParas(rmID);
          PSO_OPTS[LOOPID + 1].varUpper.set_size(OBJ.dParas(rmID));
          PSO_OPTS[LOOPID + 1].varUpper = OBJ.parasUpper.submat(rmID, 0, rmID, OBJ.dParas(rmID) - 1);
          PSO_OPTS[LOOPID + 1].varLower.set_size(OBJ.dParas(rmID));
          PSO_OPTS[LOOPID + 1].varLower = OBJ.parasLower.submat(rmID, 0, rmID, OBJ.dParas(rmID) - 1);

          PSO_Result InnerResult;
          PSO_MAIN(LOOPID + 1, PSO_OPTS, LBFGS_OPTION, OBJ, MODEL_COLLECTOR, &PSO_EXT, FALSE, FALSE, InnerResult);
          R_PARA_tmp = InnerResult.GBest;
          eff_vals(i) = InnerResult.fGBest;
        } else {
          eff_vals(i) = minDistCalc(LBFGS_OPTION, OBJ, MODEL_COLLECTOR, i, DESIGN, WT, R_PARA_tmp);
        }
        R_PARA.submat(i+1, 0, i+1, OBJ.dParas(i+1) - 1) = R_PARA_tmp;
			}
			eff_vals = eff_vals/std_vals;
			val = eff_vals.min();
			break;
		}
	}
	return val;
}

// Minimal Distance Between Two Models
double minDistCalc(const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ, model_diff_func *MODEL_COLLECTOR[], const int PAIRID,
									 const arma::mat DESIGN, const arma::rowvec WT, arma::rowvec &R_PARA_OUT)
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
  LBFGS_EVAL.WT         = WT;
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
  CONV = lbfgs(dParas, R_PARA, &fx, evaluate, NULL, &LBFGS_EVAL, &LBFGS_PAR);
  //Rprintf("CONV: %d\n", CONV);
  for (int d = 0; d < dParas; d++) { R_PARA1[d] = R_PARA[d]; }

  int LBFGS_RETRY = LBFGS_OPTION.LBFGS_RETRY; if (CONV) { LBFGS_RETRY++; }
  if ((!std::isfinite(fx)) | std::isnan(fx)) { LBFGS_RETRY++; }
  int count = 0;
  while ((count < LBFGS_RETRY)) {
    arma::rowvec SIGN_RAND = randu(1, dParas)*2.0 - 1.0;
    for (int d = 0; d < dParas; d++) { R_PARA1[d] = SIGN_RAND(d)*R_PARA1[d]; }
    CONV = lbfgs(dParas, R_PARA1, &fx1, evaluate, NULL, &LBFGS_EVAL, &LBFGS_PAR);
    //Rprintf("CONV: %d\n", CONV);
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

// Equivalence Theorem
arma::rowvec directionalDerivative(const OBJ_INFO OBJ, const arma::mat dsGrid, const arma::mat PARA_SET, const arma::rowvec alpha,
                                   model_diff_func *MODEL_COLLECTOR[])
{
  int crit_type = OBJ.crit_type;
  arma::rowvec dirDer(dsGrid.n_rows, fill::zeros);
  switch (crit_type) {
    case 0: // Fixed True
    {
      dirDer = distCalc(OBJ, dsGrid, PARA_SET, MODEL_COLLECTOR, 0);
      break;
    }
    case 1: // Max-min, Fixed True
    {
      arma::rowvec std_vals = OBJ.std_vals;
      for (int i = 0; i < OBJ.N_PAIR; i++) {
        arma::rowvec DIV = distCalc(OBJ, dsGrid, PARA_SET, MODEL_COLLECTOR, i);
        dirDer += (alpha(i)/std_vals(i))*DIV;
      }
      break;
    }
  }
  return dirDer;
}

arma::rowvec distCalc(const OBJ_INFO OBJ, const arma::mat x, const arma::mat PARA_SET,
                      model_diff_func *MODEL_COLLECTOR[], const int PAIRID)
{
  model_diff_func* func_input = MODEL_COLLECTOR[PAIRID];
  Rcpp::EvalBase *m1_func = (Rcpp::EvalBase *) func_input->M1_FUNC;
  Rcpp::EvalBase *m2_func = (Rcpp::EvalBase *) func_input->M2_FUNC;
  Rcpp::EvalBase *distFunc = (Rcpp::EvalBase *) func_input->DISTFUNC;

  arma::imat MODEL_PAIR = OBJ.MODEL_PAIR;
  int tmID = MODEL_PAIR(PAIRID, 0);
  int rmID = MODEL_PAIR(PAIRID, 1);

  arma::rowvec T_PARA = PARA_SET.submat(tmID, 0, tmID, OBJ.dParas(tmID) - 1);
  arma::rowvec R_PARA = PARA_SET.submat(rmID, 0, rmID, OBJ.dParas(rmID) - 1);

  Shield<SEXP> DESIGN_SEXP(Rcpp::wrap(x));
  Shield<SEXP> T_PARA_SEXP(Rcpp::wrap(T_PARA));
  Shield<SEXP> R_PARA_SEXP(Rcpp::wrap(R_PARA));

  Rcpp::NumericMatrix DESIGN_Rform = Rcpp::as<Rcpp::NumericMatrix>(DESIGN_SEXP);
  Rcpp::NumericVector T_PARA_Rform = Rcpp::as<Rcpp::NumericVector>(T_PARA_SEXP);
  Rcpp::NumericVector R_PARA_Rform = Rcpp::as<Rcpp::NumericVector>(R_PARA_SEXP);

  Rcpp::NumericVector eta_T_Rform((int)x.n_rows), eta_R_Rform((int)x.n_rows), DIV_Rform((int)x.n_rows);
  eta_T_Rform = (Rcpp::NumericVector) m1_func->eval(DESIGN_Rform, T_PARA_Rform);
  eta_R_Rform = (Rcpp::NumericVector) m2_func->eval(DESIGN_Rform, R_PARA_Rform);
  DIV_Rform = (Rcpp::NumericVector) distFunc->eval(eta_T_Rform, eta_R_Rform);

  arma::rowvec DIV(DIV_Rform.begin(), DIV_Rform.size());
  DIV.elem(find_nonfinite(DIV)).fill(1e10);

  return DIV;
}
