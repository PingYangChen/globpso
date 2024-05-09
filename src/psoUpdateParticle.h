// DECLARE FUNCTIONS
arma::mat expTail(arma::mat u, double b);

// BODY
void psoUpdateParticle(PSO_OPTIONS PSO_OPTS, const PSO_DYN PSO_DYN, 
                       const arma::vec fSwarm, const arma::mat PBest, const arma::rowvec GBest, 
							 				 const arma::rowvec velMax, const arma::rowvec varUpper, const arma::rowvec varLower,
							 				 arma::mat &vStep, arma::mat &swarm)
{
  int dSwarm = (int)swarm.n_cols;
  int nSwarm = (int)swarm.n_rows;
	int typePSO = PSO_OPTS.typePSO;
	double c1 = PSO_OPTS.c1;
  double c2 = PSO_OPTS.c2;
  //double chi = PSO_OPTS.chi;
	
	arma::mat velMax_Mat = repmat(velMax, nSwarm, 1);
	arma::mat GBmat = repmat(GBest, nSwarm, 1);
	GetRNGstate();
	switch (typePSO) {
		case 0: // Linearly Decreasing Weight PSO (Shi, Y. H.	and Eberhart, R. C., 1998)
		{	// The most common one
			vStep = PSO_DYN.w_cur*vStep + c1*arma::randu(nSwarm, dSwarm) % (PBest - swarm) + 
																		c2*arma::randu(nSwarm, dSwarm) % (GBmat - swarm);
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
		*/
		case 2: // Quantum PSO (Sun, J., Feng, B. and Xu, W., 2004)
		{
			arma::mat R1 = arma::randu(nSwarm, dSwarm); arma::mat R2 = arma::randu(nSwarm, dSwarm);
			arma::mat PHI = (c1*R1)/(c1*R1 + c2*R2);
			arma::mat PM = PHI % PBest + (1.0 - PHI) % GBmat;
			
			arma::mat QuantumMove;
			if (PSO_OPTS.Q_cen_type == 0) {
				QuantumMove = PSO_DYN.Q_a_cur*arma::abs(PM - swarm) % arma::log(1/arma::randu(nSwarm, dSwarm));
			} else {
				arma::mat MBest = repmat(arma::sum(PBest, 0)/((double)nSwarm), nSwarm, 1);
				QuantumMove = PSO_DYN.Q_a_cur*arma::abs(MBest - swarm) % arma::log(1/arma::randu(nSwarm, dSwarm));
			}
			
			arma::mat DICE = arma::randu(nSwarm, dSwarm);

			swarm = PM;
			swarm.elem(find(DICE > 0.5)) += QuantumMove.elem(find(DICE > 0.5));
			swarm.elem(find(DICE <= 0.5)) -= QuantumMove.elem(find(DICE <= 0.5));			
			break;
		}
		case 3: // LcRiPSO (Bonyadi, M. R., Michalewicz, Z., 2014)
		{
			arma::mat normPB = PBest + arma::randn(nSwarm, dSwarm) % repmat(PSO_DYN.LcRi_sigP, 1, dSwarm);
			arma::mat normGB = GBmat + arma::randn(nSwarm, dSwarm) % repmat(PSO_DYN.LcRi_sigG, 1, dSwarm);
			vStep = PSO_DYN.w_cur*vStep + c1*arma::randu(nSwarm, dSwarm) % (normPB - swarm) + 
																		c2*arma::randu(nSwarm, dSwarm) % (normGB - swarm);
			vStep = min(vStep, velMax_Mat); vStep = max(vStep, (-1)*velMax_Mat);
			swarm += vStep;
			break;
		}
		/*
	  case 4: // PrScPSO (Silva, A.,	Neves, A.	and Goncalves, T., 2012)
    {	// This one does not perform well.  It's needed to double-check the codes again.
      arma::mat PREY = repmat(swarm.row(PSO_DYN.PrSc_PREY), nSwarm, 1);
      arma::mat PERTURB = (repmat(varUpper, nSwarm - 1, 1) - repmat(varLower, nSwarm - 1, 1)) % (2.0*randu(nSwarm - 1, dSwarm) - 1.0);
      arma::mat PrSc_FEAR = PSO_OPTS[LOOPID].PrSc_P * arma::exp(-arma::abs(swarm.rows(1, nSwarm - 1) - repmat(swarm.row(0), nSwarm - 1, 1)));
      arma::mat DICE = randu(nSwarm - 1, dSwarm);
      PERTURB.elem(find(DICE < PrSc_FEAR)).fill(0);
      vStep.row(0) = PSO_DYN.w_cur*vStep.row(0) + 
      c1*randu(1, dSwarm) % (PREY.row(0) - swarm.row(0)) + 
      c2*randu(1, dSwarm) % (GBmat.row(0) - swarm.row(0));
      vStep.rows(1, nSwarm - 1) = PSO_DYN.w_cur*vStep.rows(1, nSwarm - 1) + 
      c1*randu(nSwarm - 1, dSwarm) % (PBest.rows(1, nSwarm - 1) - swarm.rows(1, nSwarm - 1)) + 
      c2*randu(nSwarm - 1, dSwarm) % (GBmat.rows(1, nSwarm - 1) - swarm.rows(1, nSwarm - 1)) +
      PERTURB;
      vStep = min(vStep, velMax_Mat); vStep = max(vStep, (-1)*velMax_Mat);
      swarm += vStep;
      break;
    }
		*/
  	case 5: // CSO
		{
		  //Rprintf("0\n");
		  arma::rowvec GMean = arma::sum(swarm, 0)/((double)nSwarm);
		  //Rprintf("1\n");
		  arma::uvec COUPLE_ID = arma::repmat(arma::linspace<arma::uvec>(0, nSwarm/2 - 1, nSwarm/2), 2, 1);
		  arma::vec RANDNUM = arma::randu(nSwarm, 1);
		  arma::uvec RANDIDX = sort_index(RANDNUM);
		  COUPLE_ID = COUPLE_ID(RANDIDX);
		  //Rprintf("2\n");
		  arma::uword WINNER, LOSER;
		  for (arma::uword i = 0; i < ((arma::uword)nSwarm/2); i++) {
		    //Rprintf("3\n");
		    arma::uvec COMPPAIR = find(COUPLE_ID == i);
		    if (fSwarm(COMPPAIR(0)) < fSwarm(COMPPAIR(1))) {
		      //Rprintf("4\n");
		      WINNER = COMPPAIR(0); LOSER = COMPPAIR(1);
		    } else {
		      //Rprintf("4\n");
		      WINNER = COMPPAIR(1); LOSER = COMPPAIR(0);
		    }
		    //Rprintf("5\n");
		    vStep.row(LOSER) = arma::randu(1, dSwarm) % vStep.row(LOSER) + 
		      arma::randu(1, dSwarm) % (swarm.row(WINNER) - swarm.row(LOSER)) + 
		      PSO_OPTS.CSO_phi*(arma::randu(1, dSwarm) % (GMean - swarm.row(LOSER)));
		    vStep.row(LOSER) = min(vStep.row(LOSER), velMax); 
		    vStep.row(LOSER) = max(vStep.row(LOSER), (-1)*velMax);
		    swarm.row(LOSER) += vStep.row(LOSER);
		  }
		  break;
		}
		/*
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
  	case 2024: // DExPSO
		{
		  double TE_b = PSO_OPTS.TE_b;
		  arma::mat rdPB = expTail(arma::randu(nSwarm, dSwarm), TE_b) % (PBest - swarm);
		  arma::mat rdGB = expTail(arma::randu(nSwarm, dSwarm), TE_b) % (GBmat - swarm);
		  vStep = PSO_DYN.w_cur*vStep + c1*rdPB + c2*rdGB;
		  vStep = min(vStep, velMax_Mat); vStep = max(vStep, (-1)*velMax_Mat);
		  swarm += vStep;
		  break;
		}
	  case 20241: // DExQPSO
		{
		  double TE_b = PSO_OPTS.TE_b;
		  arma::mat R1 = expTail(arma::randu(nSwarm, dSwarm), TE_b); 
		  arma::mat R2 = expTail(arma::randu(nSwarm, dSwarm), TE_b);
		  arma::mat PHI = (c1*R1)/(c1*R1 + c2*R2);
		  arma::mat PM = PHI % PBest + (1.0 - PHI) % GBmat;
		  arma::mat QuantumMove;
		  if (PSO_OPTS.Q_cen_type == 0) {
		    QuantumMove = PSO_DYN.Q_a_cur*arma::abs(PM - swarm) % arma::log(1/arma::randu(nSwarm, dSwarm));
		  } else {
		    arma::mat MBest = repmat(arma::sum(PBest, 0)/((double)nSwarm), nSwarm, 1);
		    QuantumMove = PSO_DYN.Q_a_cur*arma::abs(MBest - swarm) % arma::log(1/arma::randu(nSwarm, dSwarm));
		  }
		  arma::mat DICE = arma::randu(nSwarm, dSwarm);
		  swarm = PM;
		  swarm.elem(find(DICE > 0.5)) += QuantumMove.elem(find(DICE > 0.5));
		  swarm.elem(find(DICE <= 0.5)) -= QuantumMove.elem(find(DICE <= 0.5));       
		  break;
		}
	  case 20242: // DExCSO
		{
		  double TE_b = PSO_OPTS.TE_b;
		  arma::rowvec GMean = arma::sum(swarm, 0)/((double)nSwarm);
		  arma::uvec COUPLE_ID = arma::repmat(arma::linspace<arma::uvec>(0, nSwarm/2 - 1, nSwarm/2), 2, 1);
		  arma::vec RANDNUM = arma::randu(nSwarm, 1);
		  arma::uvec RANDIDX = sort_index(RANDNUM);
		  COUPLE_ID = COUPLE_ID(RANDIDX);
		  arma::uword WINNER, LOSER;
		  for (arma::uword i = 0; i < ((arma::uword)nSwarm/2); i++) {
		    arma::uvec COMPPAIR = find(COUPLE_ID == i);
		    if (fSwarm(COMPPAIR(0)) < fSwarm(COMPPAIR(1))) {
		      WINNER = COMPPAIR(0); LOSER = COMPPAIR(1);
		    } else {
		      WINNER = COMPPAIR(1); LOSER = COMPPAIR(0);
		    }
		    arma::mat R1 = expTail(arma::randu(1, dSwarm), TE_b);
		    arma::mat R2 = expTail(arma::randu(1, dSwarm), TE_b);
		    arma::mat R3 = expTail(arma::randu(1, dSwarm), TE_b);
		    vStep.row(LOSER) = R1 % vStep.row(LOSER) + 
		      R2 % (swarm.row(WINNER) - swarm.row(LOSER)) + 
		      PSO_OPTS.CSO_phi*(R3 % (GMean - swarm.row(LOSER)));
		    vStep.row(LOSER) = min(vStep.row(LOSER), velMax); 
		    vStep.row(LOSER) = max(vStep.row(LOSER), (-1)*velMax);
		    swarm.row(LOSER) += vStep.row(LOSER);
		  }
		  break;
		}	    
	}
	PutRNGstate();
}

arma::mat expTail(arma::mat u, double b) 
{
  double b2 = b + 2.0;
  arma::mat x = u*b2 - 1.0;
  
  arma::uvec q_left  = arma::find(u < (1.0/(b + 2.0)));
  //arma::uvec q_mid   = arma::find((u >= (1.0/(b + 2.0))) && (u <= ((b + 1.0)/(b + 2.0))));
  arma::uvec q_right = arma::find(u > ((b + 1.0)/(b + 2.0)));
  
  x.elem(q_left)  = arma::log(u.elem(q_left)*b2);
  //x.elem(q_mid)   = u.elem(q_mid)*b2 - 1.0;
  x.elem(q_right) = b - arma::log((1.0 - u.elem(q_right))*b2);
  
  /*
   if (u < (1.0/(b + 2.0))) {
   x = std::log(u*(b + 2.0));
   } else if (u > ((b + 1.0)/(b + 2.0))) {
   x = b - std::log((1.0 - u)*(b + 2.0));
   } else {
   x = u*(b + 2.0)  - 1.0;
   }
   */
  return x/b;	
}
