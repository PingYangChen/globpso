

// DEFINE STUCTURES OF PSO INFORMATION
typedef struct {
	/* Copy PSO parameters to 'genCppCode_PSOParaSetting.r'
		 for generating required C++ and R codes.
		 Please DO follow the format below.
	   type1 name1; // default1
	   type2 name2; // default2
	   type3 name3; // default3
	*/
	int seed;
	// Basic Settings
	int nSwarm; // 64
	int dSwarm; // 2
  rowvec varUpper; // c(1,1)
  rowvec varLower; // c(0,0)
  int hasInitSwarm;
  mat initSwarm;
  rowvec fixedDims;
	int maxIter; // 100
	//int checkConv; // 0
	int typePSO; // 0
  double freeRun; // 1
 	double tol; // 1e-6
 	// Basic PSO Parameters
 	double c1; // 2.05
 	double c2; // 2.05
 	double w0; // 1.2
 	double w1; // 0.2
 	double w_var; // 0.8
 	//double chi; // 0.729
 	double vk; // 4
 	// Guarantee Convergence PSO Parameters
	//int GC_S_ROOF; // 5
	//int GC_F_ROOF; // 15
	//double GC_RHO; // 1
	// Quantum PSO Parameters
  int Q_cen_type; // 0
	double Q_a0; // 1.7
	double Q_a1; // 0.7
	double Q_a_var; // 0.8
  // LcRiPSO
  double LcRi_L; // 0.01
  // CSO
  double CSO_phi; // 0.1
  // DExPSO
  double TE_b; // 2
} PSO_OPTIONS, *Ptr_PSO_OPTIONS;

// DEFINE STUCTURES OF PSO PARAMETERS WHICH WILL CHANGE ITERATIVELY
typedef struct {
  //int succ_GB; // 0
  // inertia weight
  int w_varyfor; // (int)(w_var*maxIter);
  double w_cur; // w0;
  double w_dec; // (w0 - w1)/w_varyfor;
  // GCPSO
  //int GC_S_COUNT; // 0;
  //int GC_F_COUNT; // 0;
  //double GC_RHO;
  // Quantum PSO
  int Q_a_varyfor; // (int)(Q_a_var*maxIter);
  double Q_a_cur; // Q_a0;
  double Q_a_dec; // (Q_a0 - Q_a1)/Q_a_varyfor;
  // LcRiPSO
  vec LcRi_sigP;
  vec LcRi_sigG;
} PSO_DYN, *Ptr_PSO_DYN;

// DEFINE PSO RESULTS
typedef struct {
  arma::rowvec GBest;
  double fGBest;
  arma::rowvec fGBestHist;
  arma::mat PBest;
  arma::vec fPBest;
} PSO_Result, *Ptr_PSO_Result;


// DECLARE FUNCTIONS
void psoFuncEval(const bool IF_PARALLEL, Rcpp::EvalBase *objfunc, const arma::mat swarm, arma::vec &fSwarm);
void getAlgStruct(PSO_OPTIONS &PSO_OPTS, const Rcpp::List PSO_INFO_LIST);
void PSO_MAIN(PSO_OPTIONS PSO_OPTS, Rcpp::EvalBase *objfunc,
              const bool IF_PARALLEL, const bool COUNTER_ON, PSO_Result &PSO_Result);
void psoUpdateParticle(PSO_OPTIONS PSO_OPTS, const PSO_DYN PSO_DYN, 
                       const arma::vec fSwarm, const arma::mat PBest, const arma::rowvec GBest, 
							 				 const arma::rowvec velMax, const arma::rowvec varUpper, const arma::rowvec varLower,
							 				 arma::mat &vStep, arma::mat &swarm);
void psoCheckParticle(PSO_OPTIONS PSO_OPTS, const PSO_DYN PSO_DYN, 
							 				const arma::rowvec varUpper, const arma::rowvec varLower, arma::mat &swarm);
void psoUpdateDynPara(PSO_OPTIONS PSO_OPTS, const int iter, PSO_DYN &PSO_DYN, 
											const arma::mat swarm, const arma::mat PBest, const arma::rowvec GBest,
											const arma::vec fSwarm, const arma::vec fPBest, const double fGBest);
void psoFuncEval(const bool IF_PARALLEL, Rcpp::EvalBase *objfunc, const arma::mat swarm, arma::vec &fSwarm);

#include "psoCheckParticle.h"
#include "psoUpdateParticle.h"
#include "psoUpdateDynPara.h"
#include "psoKernel.h"

// BODY

// PSO OPTIONS
void getAlgStruct(PSO_OPTIONS &PSO_OPTS, const Rcpp::List PSO_INFO_LIST)
{
  //PSO_OPTS.seed       = (int)Rcpp::as<int>(PSO_INFO_LIST["seed"]);
  PSO_OPTS.nSwarm     = (int)Rcpp::as<int>(PSO_INFO_LIST["nSwarm"]);
  PSO_OPTS.dSwarm     = (int)Rcpp::as<int>(PSO_INFO_LIST["dSwarm"]);

  Rcpp::NumericVector varUpper_Tmp   = Rcpp::as<Rcpp::NumericVector>(PSO_INFO_LIST["varUpper"]);
  arma::rowvec varUpper(varUpper_Tmp.begin(), varUpper_Tmp.size(), false);
  Rcpp::NumericVector varLower_Tmp   = Rcpp::as<Rcpp::NumericVector>(PSO_INFO_LIST["varLower"]);
  arma::rowvec varLower(varLower_Tmp.begin(), varLower_Tmp.size(), false);
  PSO_OPTS.varUpper   = varUpper;
  PSO_OPTS.varLower   = varLower;

  PSO_OPTS.hasInitSwarm = (int)Rcpp::as<int>(PSO_INFO_LIST["hasInitSwarm"]);
  if (PSO_OPTS.hasInitSwarm > 0) {
    Rcpp::NumericMatrix initSwarm_Tmp   = Rcpp::as<Rcpp::NumericMatrix>(PSO_INFO_LIST["initSwarm"]);
    arma::mat initSwarm(initSwarm_Tmp.begin(), initSwarm_Tmp.nrow(), initSwarm_Tmp.ncol(), false);
    PSO_OPTS.initSwarm  = initSwarm;
  }

  Rcpp::NumericVector fixedDims_Tmp = Rcpp::as<Rcpp::NumericVector>(PSO_INFO_LIST["fixedDims"]);
  arma::rowvec fixedDims(fixedDims_Tmp.begin(), fixedDims_Tmp.size(), false);
  PSO_OPTS.fixedDims = fixedDims;
  
  PSO_OPTS.maxIter    = (int)Rcpp::as<int>(PSO_INFO_LIST["maxIter"]);
  //PSO_OPTS.checkConv  = (int)Rcpp::as<int>(PSO_INFO_LIST["checkConv"]);
  PSO_OPTS.typePSO    = (int)Rcpp::as<int>(PSO_INFO_LIST["typePSO"]);
  PSO_OPTS.freeRun    = (double)Rcpp::as<double>(PSO_INFO_LIST["freeRun"]);
  PSO_OPTS.tol        = (double)Rcpp::as<double>(PSO_INFO_LIST["tol"]);
  PSO_OPTS.c1         = (double)Rcpp::as<double>(PSO_INFO_LIST["c1"]);
  PSO_OPTS.c2         = (double)Rcpp::as<double>(PSO_INFO_LIST["c2"]);
  PSO_OPTS.w0         = (double)Rcpp::as<double>(PSO_INFO_LIST["w0"]);
  PSO_OPTS.w1         = (double)Rcpp::as<double>(PSO_INFO_LIST["w1"]);
  PSO_OPTS.w_var      = (double)Rcpp::as<double>(PSO_INFO_LIST["w_var"]);
  //PSO_OPTS.chi        = (double)Rcpp::as<double>(PSO_INFO_LIST["chi"]);
  PSO_OPTS.vk         = (double)Rcpp::as<double>(PSO_INFO_LIST["vk"]);
  //PSO_OPTS.GC_S_ROOF  = (int)Rcpp::as<int>(PSO_INFO_LIST["GC_S_ROOF"]);
  //PSO_OPTS.GC_F_ROOF  = (int)Rcpp::as<int>(PSO_INFO_LIST["GC_F_ROOF"]);
  //PSO_OPTS.GC_RHO     = (double)Rcpp::as<double>(PSO_INFO_LIST["GC_RHO"]);
  PSO_OPTS.Q_cen_type = (int)Rcpp::as<int>(PSO_INFO_LIST["Q_cen_type"]);
  PSO_OPTS.Q_a0       = (double)Rcpp::as<double>(PSO_INFO_LIST["Q_a0"]);
  PSO_OPTS.Q_a1       = (double)Rcpp::as<double>(PSO_INFO_LIST["Q_a1"]);
  PSO_OPTS.Q_a_var    = (double)Rcpp::as<double>(PSO_INFO_LIST["Q_a_var"]);
  PSO_OPTS.LcRi_L     = (double)Rcpp::as<double>(PSO_INFO_LIST["LcRi_L"]);
  PSO_OPTS.CSO_phi    = (double)Rcpp::as<double>(PSO_INFO_LIST["CSO_phi"]);
  PSO_OPTS.TE_b    = (double)Rcpp::as<double>(PSO_INFO_LIST["TE_b"]);
}

