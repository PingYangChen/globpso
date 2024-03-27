// DECLARE FUNCTIONS

// BODY
void deUpdatePop(DE_OPTIONS DE_OPTS, const DE_DYN DE_DYN, 
                 const arma::vec fpop, const arma::mat PBest, const arma::rowvec GBest, 
							 	 const arma::rowvec varUpper, const arma::rowvec varLower,
							 	 const arma::mat &pop, arma::mat &vMuta, arma::mat &newpop)
{
  int dPop = (int)pop.n_cols;
  int nPop = (int)pop.n_rows;
	int typeDE = DE_OPTS.typeDE;
	double sf = DE_OPTS.sf;
  double cr = DE_OPTS.cr;

	// arma::mat GBmat = repmat(GBest, npop, 1);
	switch (typeDE) {
		case 0: // Default Differential Eolution (Storn, R. and Price, K., 1997)
		{	// The most common one
      arma::umat rand_pop(3, nPop);
      arma::uvec rand_tmp(3);
      for (int i = 0; i < nPop; i++) {
        rand_tmp = arma::randperm(nPop - 1, 3);
        rand_tmp.elem(arma::find(rand_tmp >= i)) += 1;
        rand_pop.col(i) = rand_tmp;
      }
      vMuta = pop.rows(rand_pop.row(0)) + sf*(pop.rows(rand_pop.row(1)) - pop.rows(rand_pop.row(2)));
      for (int i = 0; i < nPop; i++) {
        int rand_c = arma::as_scalar(arma::randi(1, distr_param(0, dPop - 1)));
        for (int j = 0; j < dPop; j++) {
          double coin = arma::as_scalar(arma::randu(1));
          if ((coin <= cr) | (j == rand_c)) {
            newpop(i, j) = vMuta(i, j);
          } else {
            newpop(i, j) = pop(i, j);
          }
        }
      }
			break;
		}
	}	
}
