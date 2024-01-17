void psoUpdateDynPara(PSO_OPTIONS PSO_OPTS, const int iter, PSO_DYN &PSO_DYN, 
											const arma::mat swarm, const arma::mat PBest, const arma::rowvec GBest,
											const arma::vec fSwarm, const arma::vec fPBest, const double fGBest)
{
	int typePSO = PSO_OPTS.typePSO;
  if (iter < 0) { // INITIALIZE
    // The most common one
    // Linearly Decreasing Weight PSO (Shi, Y. H.	and Eberhart, R. C., 1998)
    int w_varyfor = (int)(PSO_OPTS.w_var*PSO_OPTS.maxIter);
    PSO_DYN.w_varyfor	= w_varyfor;
    PSO_DYN.w_cur			= PSO_OPTS.w0;
    PSO_DYN.w_dec			= (PSO_OPTS.w0 - PSO_OPTS.w1)/w_varyfor;      // Inertia weight change per iteration step
		if (typePSO == 2) { // Quantum PSO (Sun, J., Feng, B. and Xu, W., 2004)
  		int Q_a_varyfor = (int)(PSO_OPTS.Q_a_var*PSO_OPTS.maxIter);
	  	PSO_DYN.Q_a_varyfor	= Q_a_varyfor;
			PSO_DYN.Q_a_cur			= PSO_OPTS.Q_a0;
			PSO_DYN.Q_a_dec			= (PSO_OPTS.Q_a0 - PSO_OPTS.Q_a1)/Q_a_varyfor;
    }
		if (typePSO == 4) { // LcRiPSO
	    arma::mat GBmat = repmat(GBest, PSO_OPTS.nSwarm, 1);
	    arma::vec EUD_PB = arma::sqrt(sum((swarm - PBest) % (swarm - PBest), 1));
	    arma::vec EUD_GB = arma::sqrt(sum((swarm - GBmat) % (swarm - GBmat), 1));
	    arma::vec LcRi_sigP(PSO_OPTS.nSwarm); LcRi_sigP.fill(EUD_PB.max());
	    arma::vec LcRi_sigG(PSO_OPTS.nSwarm); LcRi_sigG.fill(EUD_GB.max());
	    PSO_DYN.LcRi_sigP = LcRi_sigP;
	    PSO_DYN.LcRi_sigG = LcRi_sigG;
	  }
		/*
		// Guarantee Convergence PSO 
		//PSO_DYN.succ_GB	= 0;
		PSO_DYN.GC_S_COUNT	= 0; 
		PSO_DYN.GC_F_COUNT	= 0;
		PSO_DYN.GC_RHO			= PSO_OPTS.GC_RHO; 
    // PrScPSO
    arma::vec PREY_SET = fSwarm.subvec(1, PSO_OPTS[LOOPID].nSwarm - 1);
    arma::uword PREY_ID; 
    PREY_SET.min(PREY_ID); 
    PSO_DYN.PrSc_PREY = PREY_ID + 1;
    PSO_DYN.PrSc_SIG = EUD_GB.max();
		*/
  } else { // UPDATE
		// Linearly Decreasing Weight PSO (Shi, Y. H.	and Eberhart, R. C., 1998)
		if (iter <= PSO_DYN.w_varyfor) {
		  PSO_DYN.w_cur = PSO_DYN.w_cur - PSO_DYN.w_dec;
		}
		if (typePSO == 2) {// Quantum PSO (Sun, J., Feng, B. and Xu, W., 2004)
    	if (iter <= PSO_DYN.Q_a_varyfor) {
    	  PSO_DYN.Q_a_cur = PSO_DYN.Q_a_cur - PSO_DYN.Q_a_dec;
    	}
		}
  	if (typePSO == 4) { // LcRiPSO
  	  arma::mat GBmat = repmat(GBest, PSO_OPTS.nSwarm, 1);
  	  arma::vec EUD_PB = arma::sqrt(sum((swarm - PBest) % (swarm - PBest), 1));
  	  arma::vec EUD_GB = arma::sqrt(sum((swarm - GBmat) % (swarm - GBmat), 1));
  	  arma::uvec EUD_ZERO_PB = find(EUD_PB == 0);
  	  arma::uvec EUD_ZERO_GB = find(EUD_GB == 0);
  	  EUD_PB.elem(EUD_ZERO_PB) = PSO_DYN.LcRi_sigP(EUD_ZERO_PB);
  	  EUD_GB.elem(EUD_ZERO_GB) = PSO_DYN.LcRi_sigG(EUD_ZERO_GB);
  	  PSO_DYN.LcRi_sigP = PSO_OPTS.LcRi_L * EUD_PB;
  	  PSO_DYN.LcRi_sigG = PSO_OPTS.LcRi_L * EUD_GB;
  	}
    /*
		// GCPSO
		if (PSO_OPTS.typePSO == 1) {
			if (PSO_DYN.succ_GB == 1) { 
				PSO_DYN.GC_S_COUNT++; PSO_DYN.GC_F_COUNT = 0; PSO_DYN.succ_GB = 0;
			} else {
				PSO_DYN.GC_F_COUNT++; PSO_DYN.GC_S_COUNT = 0; 
			}
			if (PSO_DYN.GC_S_COUNT > PSO_OPTS.GC_S_ROOF) {
				PSO_DYN.GC_RHO *= 2.0; //Rprintf("wider\n");
			} else if (PSO_DYN.GC_F_COUNT > PSO_OPTS.GC_F_ROOF) {
				PSO_DYN.GC_RHO *= 0.5; //Rprintf("narrower\n");
			}
		}
    // PrScPSO
    arma::vec PREY_SET = fSwarm.subvec(1, PSO_OPTS[LOOPID].nSwarm - 1);
    arma::uword PREY_ID; 
    PREY_SET.min(PREY_ID); 
    PSO_DYN.PrSc_PREY = PREY_ID + 1;
    PSO_DYN.PrSc_SIG = EUD_GB.max();
    */
  }
}
