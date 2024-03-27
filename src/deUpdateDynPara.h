void deUpdateDynPara(DE_OPTIONS DE_OPTS, const int iter, DE_DYN &DE_DYN, 
										 const arma::mat pop, const arma::mat PBest, const arma::rowvec GBest,
										 const arma::vec fpop, const arma::vec fPBest, const double fGBest)
{
	//int typeDE = DE_OPTS.typeDE;
  if (iter < 0) { // INITIALIZE
    // The most common one
   
  } else { // UPDATE
		//
		/*
		if (iter <= DE_DYN.sf_varyfor) {
		  DE_DYN.sf_cur = DE_DYN.sf_cur - DE_DYN.sf_dec;
		}
		*/
		
  }
}
