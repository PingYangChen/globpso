void psoUpdateDynPara(const int LOOPID, PSO_OPTIONS PSO_OPTS[], const int iter, PSO_DYN &PSO_DYN, 
											const arma::mat swarm, const arma::mat PBest, const arma::rowvec GBest,
											const arma::vec fSwarm, const arma::vec fPBest, const double fGBest)
{
  if (iter < 0) { // INITIALIZE

  	//PSO_DYN.succ_GB	= 0;

  	int w_varyfor = (int)(PSO_OPTS[LOOPID].w_var*PSO_OPTS[LOOPID].maxIter);
  	PSO_DYN.w_varyfor	= w_varyfor;
	  PSO_DYN.w_cur			= PSO_OPTS[LOOPID].w0;
	  PSO_DYN.w_dec			= (PSO_OPTS[LOOPID].w0 - PSO_OPTS[LOOPID].w1)/w_varyfor;      // Inertia weight change per iteration step
	  // Quantum PSO 
	  /*
	  int Q_a_varyfor = (int)(PSO_OPTS[LOOPID].Q_a_var*PSO_OPTS[LOOPID].maxIter);
  	PSO_DYN.Q_a_varyfor	= Q_a_varyfor;
		PSO_DYN.Q_a_cur			= PSO_OPTS[LOOPID].Q_a0;
		PSO_DYN.Q_a_dec			= (PSO_OPTS[LOOPID].Q_a0 - PSO_OPTS[LOOPID].Q_a1)/Q_a_varyfor;
		// Guarantee Convergence PSO 
		PSO_DYN.GC_S_COUNT	= 0; 
		PSO_DYN.GC_F_COUNT	= 0;
		PSO_DYN.GC_RHO			= PSO_OPTS[LOOPID].GC_RHO; 
		// LcRiPSO
		arma::mat GBmat = repmat(GBest, PSO_OPTS[LOOPID].nSwarm, 1);
		arma::vec EUD_PB = arma::sqrt(sum((swarm - PBest) % (swarm - PBest), 1));
		arma::vec EUD_GB = arma::sqrt(sum((swarm - GBmat) % (swarm - GBmat), 1));
		arma::vec LcRi_sigP(PSO_OPTS[LOOPID].nSwarm); LcRi_sigP.fill(EUD_PB.max());
		arma::vec LcRi_sigG(PSO_OPTS[LOOPID].nSwarm); LcRi_sigG.fill(EUD_GB.max());
		PSO_DYN.LcRi_sigP = LcRi_sigP;
		PSO_DYN.LcRi_sigG = LcRi_sigG;
		*/
  } else { // UPDATE
    
		if (iter <= PSO_DYN.w_varyfor) 		PSO_DYN.w_cur = PSO_DYN.w_cur - PSO_DYN.w_dec; 
		/*
		// QPSO
    if (iter <= PSO_DYN.Q_a_varyfor)	PSO_DYN.Q_a_cur = PSO_DYN.Q_a_cur - PSO_DYN.Q_a_dec; 
		// GCPSO
		if (PSO_OPTS[LOOPID].typePSO == 1) {
			if (PSO_DYN.succ_GB == 1) { 
				PSO_DYN.GC_S_COUNT++; PSO_DYN.GC_F_COUNT = 0; PSO_DYN.succ_GB = 0;
			} else {
				PSO_DYN.GC_F_COUNT++; PSO_DYN.GC_S_COUNT = 0; 
			}
			if (PSO_DYN.GC_S_COUNT > PSO_OPTS[LOOPID].GC_S_ROOF) {
				PSO_DYN.GC_RHO *= 2.0; //Rprintf("wider\n");
			} else if (PSO_DYN.GC_F_COUNT > PSO_OPTS[LOOPID].GC_F_ROOF) {
				PSO_DYN.GC_RHO *= 0.5; //Rprintf("narrower\n");
			}
		}
		// LcRiPSO
		arma::mat GBmat = repmat(GBest, PSO_OPTS[LOOPID].nSwarm, 1);
		arma::vec EUD_PB = arma::sqrt(sum((swarm - PBest) % (swarm - PBest), 1));
		arma::vec EUD_GB = arma::sqrt(sum((swarm - GBmat) % (swarm - GBmat), 1));
		arma::uvec EUD_ZERO_PB = find(EUD_PB == 0);
		arma::uvec EUD_ZERO_GB = find(EUD_GB == 0);
		EUD_PB.elem(EUD_ZERO_PB) = PSO_DYN.LcRi_sigP(EUD_ZERO_PB);
		EUD_GB.elem(EUD_ZERO_GB) = PSO_DYN.LcRi_sigG(EUD_ZERO_GB);
		PSO_DYN.LcRi_sigP = PSO_OPTS[LOOPID].LcRi_L * EUD_PB;
		PSO_DYN.LcRi_sigG = PSO_OPTS[LOOPID].LcRi_L * EUD_GB;
		*/
  }
}

/* Under Construction Now
arma::mat getTopology(const int LOOPID, const PSO_OPTIONS PSO_OPTS[], 
						 					const arma::mat swarm, const arma::mat PBest, const arma::rowvec GBest,
											const arma::vec fSwarm, const arma::vec fPBest, const double fGBest)
{
	int typeTopo = PSO_OPTS[LOOPID].typeTopo;
	int nSwarm = PSO_OPTS[LOOPID].nSwarm;
	int nGroup = PSO_OPTS[LOOPID].nGroup;
	arma::mat TOPOMAT;
	switch (typeTopo) {
		case 0: { // ALL
			TOPOMAT.set_size(nSwarm, 1); TOPOMAT.fill(1);
			break;
		}
		case 1: { // RAND
			TOPOMAT.set_size(nSwarm, nGroup);
			arma::vec RANDVAL = randu(nSwarm, 1);
			arma::uvec RANDIDX = sort_index(RANDVAL);
			int minMember = (int)(std::round(((double)nSwarm)/(double)nGroup));
			arma::vec TOPOTMP(nSwarm);
			for (int i = 0; i < nGroup; i++) {
				TOPOTMP.zeros();
				TOPOTMP(RANDIDX.subvec(i*minMember, (i+1)*minMember - 1)).fill(1);
				TOPOMAT.col(i) = TOPOTMP;
			}
			break;
		}
		case 2: { //RING
			TOPOMAT.set_size(nSwarm, nSwarm); TOPOMAT.zeros();
			arma::vec TOPOTMP(nSwarm);
			for (int i = 0; i < nSwarm; i++) {
				TOPOTMP.zeros();
				if (i == 0) { TOPOTMP(nSwarm - 1) = 1; TOPOTMP.subvec(0,1).fill(1); }
				else if (i == (nSwarm - 1)) { TOPOTMP(0) = 1; TOPOTMP.subvec(nSwarm - 2,nSwarm - 1).fill(1); }
				else TOPOTMP.subvec(i-1,i+1).fill(1);
				TOPOMAT.col(i) = TOPOTMP;
			}
			break;
		}
	}
	return TOPOMAT;
}*/