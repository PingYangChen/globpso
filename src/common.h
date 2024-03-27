// Rcpp Header File
#include <cmath>
#include <RcppArmadillo.h>
//#include <omp.h>

//using namespace Rcpp;
using namespace arma;

#ifndef Rcpp_EVALUATE_H_
#define Rcpp_EVALUATE_H_

namespace Rcpp {

  class EvalBase {
    public:
        EvalBase() : neval(0) {};
        virtual double eval(SEXP x) = 0;
        virtual ~EvalBase() {};
        //unsigned long getNbEvals() { return neval; }
    protected:
        unsigned long int neval;
  };

  class EvalStandard : public EvalBase {
    public:
        EvalStandard(SEXP fcall_, SEXP env_) : fcall(fcall_), env(env_) {}
        double eval(SEXP x) {
          //neval++;
          return defaultfun(x);
        }
        ~EvalStandard() {};
    private:
        SEXP fcall, env;
        double defaultfun(SEXP x) {
          //Shield<SEXP> fn(Rcpp::Rcpp_lang3(fcall, x, R_DotsSymbol));
          Shield<SEXP> fn(::Rf_lang3(fcall, x, R_DotsSymbol));
          Shield<SEXP> sexp_fvec(::Rf_eval(fn, env));
          //SEXP sexp_fvec = Rcpp::Rcpp_eval(fn, env); // too slow
          double f_result = (double)Rcpp::as<double>(sexp_fvec);
          return f_result;
        }
    };

  typedef double (*funcPtr)(SEXP, SEXP);

  class EvalCompiled : public EvalBase {
    public:
        EvalCompiled(Rcpp::XPtr<funcPtr> xptr, SEXP __env) {
          funptr = *(xptr);
          env = __env;
        };
        EvalCompiled(SEXP xps, SEXP __env) {
          Rcpp::XPtr<funcPtr> xptr(xps);
          funptr = *(xptr);
          env = __env;
        };
        double eval(SEXP x) {
          //neval++;
          double f_result = funptr(x, env);
          return f_result;
        }
        ~EvalCompiled() {};
    private:
        funcPtr funptr;
        SEXP env;
    };
}

#endif

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

// DECLARE FUNCTIONS

// BODY
void matrixPrintf(const mat &m)
{
  for (uword i = 0; i < m.n_rows; i++) {
    for (uword j = 0; j < m.n_cols; j++) Rprintf("%4.4f\t", m(i,j));
    Rprintf("\n");
  }
  Rprintf("\n\n");
}

void umatrixPrintf(const umat &m)
{
  for (uword i = 0; i < m.n_rows; i++) {
    for (uword j = 0; j < m.n_cols; j++) Rprintf("%d\t", m(i,j));
    Rprintf("\n");
  }
  Rprintf("\n\n");
}

void rvecPrintf(const rowvec &v)
{
  for (uword i = 0; i < v.n_elem; i++) Rprintf("%4.4f\t", v(i));
  Rprintf("\n\n");
}

void vecPrintf(const vec &v)
{
  for (uword i = 0; i < v.n_elem; i++) Rprintf("%4.4f\t", v(i));
  Rprintf("\n\n");
}

void uvecPrintf(const uvec &v)
{
  for (uword i = 0; i < v.n_elem; i++) Rprintf("%d\t", v(i));
  Rprintf("\n\n");
}
