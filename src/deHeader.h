

// DEFINE STUCTURES OF DE INFORMATION
typedef struct {
	/* Copy PSO parameters to 'genCppCode_DEParaSetting.r'
		 for generating required C++ and R codes.
		 Please DO follow the format below.
	   type1 name1; // default1
	   type2 name2; // default2
	   type3 name3; // default3
	*/
	int seed;
	// Basic Settings
	int nPop; // 64
	int dPop; // 2
  rowvec varUpper; // c(1,1)
  rowvec varLower; // c(0,0)
  int hasInitPop;
  mat initPop;
  rowvec fixedDims;
	int maxIter; // 100
	//int checkConv; // 0
	int typeDE; // 0
  double freeRun; // 1
 	double tol; // 1e-6
 	// Basic DE Parameters
 	double sf; // 0.5 scaling factor
 	double cr; // 0.1 crossover rate
} DE_OPTIONS, *Ptr_DE_OPTIONS;

// DEFINE STUCTURES OF DE PARAMETERS WHICH WILL CHANGE ITERATIVELY
typedef struct {
  // scaling factor
  int sf_varfor; // (int)(sf_var*maxIter);
  double sf_cur; //
  double sf_dec; //
  // crossover rate
  int cr_varfor; // (int)(cr_var*maxIter);
  double cr_cur; //
  double cr_dec; // (p0 - p1)/p_varyfor;
} DE_DYN, *Ptr_DE_DYN;

// DEFINE DE RESULTS
typedef struct {
  arma::rowvec GBest;
  double fGBest;
  arma::rowvec fGBestHist;
  arma::mat PBest;
  arma::vec fPBest;
} DE_Result, *Ptr_DE_Result;


// DECLARE FUNCTIONS
void getParamDE(DE_OPTIONS &DE_OPTS, const Rcpp::List DE_INFO_LIST);
void DE_MAIN(DE_OPTIONS DE_OPTS, Rcpp::EvalBase *objfunc,
             const bool IF_PARALLEL, const bool COUNTER_ON, DE_Result &DE_Result);
void deUpdatePop(DE_OPTIONS DE_OPTS, const DE_DYN DE_DYN, 
                 const arma::vec fpop, const arma::mat PBest, const arma::rowvec GBest, 
                 const arma::rowvec varUpper, const arma::rowvec varLower,
                 const arma::mat &pop, arma::mat &vMuta, arma::mat &newpop);
void deCheckPop(DE_OPTIONS DE_OPTS, const DE_DYN DE_DYN, 
                const arma::rowvec varUpper, const arma::rowvec varLower, arma::mat &pop);
void deUpdateDynPara(DE_OPTIONS DE_OPTS, const int iter, DE_DYN &DE_DYN, 
                     const arma::mat pop, const arma::mat PBest, const arma::rowvec GBest,
                     const arma::vec fpop, const arma::vec fPBest, const double fGBest);

#include "deCheckPop.h"
#include "deUpdatePop.h"
#include "deUpdateDynPara.h"
#include "deKernel.h"


// PSO OPTIONS
void getParamDE(DE_OPTIONS &DE_OPTS, const Rcpp::List DE_INFO_LIST)
{
  //DE_OPTS.seed       = (int)Rcpp::as<int>(DE_INFO_LIST["seed"]);
  DE_OPTS.nPop     = (int)Rcpp::as<int>(DE_INFO_LIST["nPop"]);
  DE_OPTS.dPop     = (int)Rcpp::as<int>(DE_INFO_LIST["dPop"]);

  Rcpp::NumericVector varUpper_Tmp   = Rcpp::as<Rcpp::NumericVector>(DE_INFO_LIST["varUpper"]);
  arma::rowvec varUpper(varUpper_Tmp.begin(), varUpper_Tmp.size(), false);
  Rcpp::NumericVector varLower_Tmp   = Rcpp::as<Rcpp::NumericVector>(DE_INFO_LIST["varLower"]);
  arma::rowvec varLower(varLower_Tmp.begin(), varLower_Tmp.size(), false);
  DE_OPTS.varUpper   = varUpper;
  DE_OPTS.varLower   = varLower;

  DE_OPTS.hasInitPop = (int)Rcpp::as<int>(DE_INFO_LIST["hasInitPop"]);
  if (DE_OPTS.hasInitPop > 0) {
    Rcpp::NumericMatrix initPop_Tmp   = Rcpp::as<Rcpp::NumericMatrix>(DE_INFO_LIST["initPop"]);
    arma::mat initPop(initPop_Tmp.begin(), initPop_Tmp.nrow(), initPop_Tmp.ncol(), false);
    DE_OPTS.initPop  = initPop;
  }

  Rcpp::NumericVector fixedDims_Tmp = Rcpp::as<Rcpp::NumericVector>(DE_INFO_LIST["fixedDims"]);
  arma::rowvec fixedDims(fixedDims_Tmp.begin(), fixedDims_Tmp.size(), false);
  DE_OPTS.fixedDims = fixedDims;
  
  DE_OPTS.maxIter    = (int)Rcpp::as<int>(DE_INFO_LIST["maxIter"]);
  DE_OPTS.typeDE     = (int)Rcpp::as<int>(DE_INFO_LIST["typeDE"]);
  DE_OPTS.freeRun    = (double)Rcpp::as<double>(DE_INFO_LIST["freeRun"]);
  DE_OPTS.tol        = (double)Rcpp::as<double>(DE_INFO_LIST["tol"]);
  DE_OPTS.sf         = (double)Rcpp::as<double>(DE_INFO_LIST["sf"]);
  DE_OPTS.cr         = (double)Rcpp::as<double>(DE_INFO_LIST["cr"]);
}


