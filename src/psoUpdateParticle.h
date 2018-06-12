
void psoUpdateParticle(const int LOOPID, PSO_OPTIONS PSO_OPTS[], const PSO_DYN PSO_DYN, 
							 				 const arma::mat PBest, const arma::rowvec GBest, 
							 				 const arma::rowvec velMax, const arma::rowvec varUpper, const arma::rowvec varLower,
							 				 arma::mat &vStep, arma::mat &swarm)
{
  int dSwarm = (int)swarm.n_cols;
  int nSwarm = (int)swarm.n_rows;
	//int typePSO = PSO_OPTS[LOOPID].typePSO;
	int typePSO = 0;
	double c1 = PSO_OPTS[LOOPID].c1;
  double c2 = PSO_OPTS[LOOPID].c2;
  //double chi = PSO_OPTS[LOOPID].chi;
	
	arma::mat velMax_Mat = repmat(velMax, nSwarm, 1);
	arma::mat GBmat = repmat(GBest, nSwarm, 1);
	switch (typePSO) {
		case 0: // Linearly Decreasing Weight PSO (Shi, Y. H.	and Eberhart, R. C., 1998)
		{	// The most common one
			vStep = PSO_DYN.w_cur*vStep + c1*randu(nSwarm, dSwarm) % (PBest - swarm) + 
																		c2*randu(nSwarm, dSwarm) % (GBmat - swarm);
			vStep = min(vStep, velMax_Mat); vStep = max(vStep, (-1)*velMax_Mat);
			swarm += vStep;													
			break;
		}
		/*
		case 1: // GCPSO (van den Bergh, F. and	Engelbrecht, A. P., 2002)
		{
			vStep.rows(0, nSwarm - 2) = PSO_DYN.w_cur*vStep.rows(0, nSwarm - 2) + 
																	c1*randu(nSwarm - 1, dSwarm) % (PBest.rows(0, nSwarm - 2) - swarm.rows(0, nSwarm - 2)) + 
																	c2*randu(nSwarm - 1, dSwarm) % (GBmat.rows(0, nSwarm - 2) - swarm.rows(0, nSwarm - 2));

			vStep.row(nSwarm - 1) = PSO_DYN.w_cur*vStep.row(nSwarm - 1) - swarm.row(nSwarm - 1) + GBest + 
															PSO_DYN.GC_RHO*(1 - 2*randu(1, dSwarm));

			vStep = min(vStep, velMax_Mat); vStep = max(vStep, (-1)*velMax_Mat);
			swarm += vStep;		
			break;
		}
		case 2: // Quantum PSO (Sun, J., Feng, B. and Xu, W., 2004)
		{
			arma::mat R1 = randu(nSwarm, dSwarm); arma::mat R2 = randu(nSwarm, dSwarm);
			arma::mat PHI = (c1*R1)/(c1*R1 + c2*R2);
			arma::mat PM = PHI % PBest + (1.0 - PHI) % GBmat;
			
			arma::mat QuantumMove;
			if (PSO_OPTS[LOOPID].Q_cen_type == 0) {
				QuantumMove = PSO_DYN.Q_a_cur*arma::abs(PM - swarm) % arma::log(1/randu(nSwarm, dSwarm));
			} else {
				arma::mat MBest = repmat(arma::sum(PBest, 0)/((double)nSwarm), nSwarm, 1);
				QuantumMove = PSO_DYN.Q_a_cur*arma::abs(MBest - swarm) % arma::log(1/randu(nSwarm, dSwarm));
			}
			
			arma::mat DICE = randu(nSwarm, dSwarm);

			swarm = PM;
			swarm.elem(find(DICE > 0.5)) += QuantumMove.elem(find(DICE > 0.5));
			swarm.elem(find(DICE <= 0.5)) -= QuantumMove.elem(find(DICE <= 0.5));			
			break;
		}
		case 3: // LcRiPSO (Bonyadi, M. R., Michalewicz, Z., 2014)
		{
			arma::mat normPB = PBest + randn(nSwarm, dSwarm) % repmat(PSO_DYN.LcRi_sigP, 1, dSwarm);
			arma::mat normGB = GBmat + randn(nSwarm, dSwarm) % repmat(PSO_DYN.LcRi_sigG, 1, dSwarm);
			vStep = PSO_DYN.w_cur*vStep + c1*randu(nSwarm, dSwarm) % (normPB - swarm) + 
																		c2*randu(nSwarm, dSwarm) % (normGB - swarm);
			vStep = min(vStep, velMax_Mat); vStep = max(vStep, (-1)*velMax_Mat);
			swarm += vStep;
			break;
		}
		case 6: // FIPSO (Mendes, R., Kennedy, J. and Neves, J., 2004)
		{

			arma::mat R1 = randu(nSwarm, dSwarm); arma::mat R2 = randu(nSwarm, dSwarm);
			arma::mat PHI = R1/(R1 + R2);
			arma::mat PM = PHI % PBest + (1.0 - PHI) % GBmat;

			arma::mat PM = repmat(arma::sum(PBest, 0)/((double)nSwarm), nSwarm, 1);

			vStep = chi*(vStep + (R1 + R2) % (PM - swarm));
			vStep = min(vStep, velMax_Mat); vStep = max(vStep, (-1)*velMax_Mat);
			swarm += vStep;

			break;
		}
		*/
	}	
}

