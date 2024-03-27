
void deCheckPop(DE_OPTIONS DE_OPTS, const DE_DYN DE_DYN, 
							 	const arma::rowvec varUpper, const arma::rowvec varLower, arma::mat &pop)
{
  int dSwarm = (int)pop.n_cols;
  int nSwarm = (int)pop.n_rows;

	arma::mat varUB_Mat = repmat(varUpper, nSwarm, 1);
	arma::mat varLB_Mat = repmat(varLower, nSwarm, 1);

	arma::mat popTmp1;
  // UPDATE POSITION
	
	arma::umat popChange;
	popChange = find(pop > varUB_Mat);
	popTmp1 = arma::pow(arma::randu(nSwarm, dSwarm), 5e-2) % (varUB_Mat - varLB_Mat) + varLB_Mat;
	for (arma::uword i = 0; i < popTmp1.n_rows; i++) {
	  for (arma::uword j = 0; j < popTmp1.n_cols; j++) {
  	  double coin = arma::as_scalar(arma::randu(1));
  	  if (coin > 0.5) {
  	    popTmp1(i, j) = varUB_Mat(i,j);
  	  }
	  }
	}
	pop.elem(popChange) = popTmp1.elem(popChange); 
			
	popChange = find(pop < varLB_Mat);
	popTmp1 = (1.0 - arma::pow(1.0 - arma::randu(nSwarm, dSwarm), 5e-2)) % (varUB_Mat - varLB_Mat) + varLB_Mat;
	for (arma::uword i = 0; i < popTmp1.n_rows; i++) {
	  for (arma::uword j = 0; j < popTmp1.n_cols; j++) {
	    double coin = arma::as_scalar(arma::randu(1));
	    if (coin > 0.5) {
	      popTmp1(i, j) = varLB_Mat(i,j);
	    }
	  }
	}
	pop.elem(popChange) = popTmp1.elem(popChange);

}
