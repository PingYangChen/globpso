
void psoCheckParticle(PSO_OPTIONS PSO_OPTS, const PSO_DYN PSO_DYN, 
							 				const arma::rowvec varUpper, const arma::rowvec varLower, arma::mat &swarm)
{
  int dSwarm = (int)swarm.n_cols;
  int nSwarm = (int)swarm.n_rows;
  //mexPrintf("P LOOP: %d\n", LOOPID);

	arma::mat varUB_Mat = repmat(varUpper, nSwarm, 1);
	arma::mat varLB_Mat = repmat(varLower, nSwarm, 1);

	arma::mat swarmTmp1;
  // UPDATE POSITION
	
	arma::umat SwarmChange;
	SwarmChange = find(swarm > varUB_Mat);
	swarmTmp1 = arma::pow(arma::randu(nSwarm, dSwarm), 1e-8) % (varUB_Mat - varLB_Mat) + varLB_Mat;
	swarm.elem(SwarmChange) = swarmTmp1.elem(SwarmChange);
			
	SwarmChange = find(swarm < varLB_Mat);
	swarmTmp1 = (1.0 - arma::pow(1.0 - arma::randu(nSwarm, dSwarm), 1e-8)) % (varUB_Mat - varLB_Mat) + varLB_Mat;
	swarm.elem(SwarmChange) = swarmTmp1.elem(SwarmChange);
	/*
	swarmTmp1 = swarm;
	swarmTmp1 = min(swarmTmp1, varUB_Mat);
	swarmTmp1 = max(swarmTmp1, varLB_Mat);
	swarm = swarmTmp1;
	 */
}
