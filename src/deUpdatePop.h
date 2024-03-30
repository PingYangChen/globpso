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
	arma::mat GBmat = repmat(GBest, nPop, 1);
	//
	GetRNGstate();
	switch (typeDE) {
		case 0: // Default Differential Evolution (DE/rand/1) (Storn, R. and Price, K., 1997)
		{	
  	  arma::uword n_rand = 3;
  	  arma::umat rand_pop(n_rand, nPop);
  	  arma::uvec rand_tmp(n_rand);
      for (int i = 0; i < nPop; i++) {
        rand_tmp = arma::randperm(nPop - 1, n_rand);
        rand_tmp.elem(arma::find(rand_tmp >= i)) += 1;
        rand_pop.col(i) = rand_tmp;
      }
      // the update formula of DE/rand/1
      vMuta = pop.rows(rand_pop.row(0)) + sf*(pop.rows(rand_pop.row(1)) - pop.rows(rand_pop.row(2)));
      // Generate new population
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
	  case 1: // DE/rand/2
		{
		  arma::uword n_rand = 5;
		  arma::umat rand_pop(n_rand, nPop);
		  arma::uvec rand_tmp(n_rand);
		  for (int i = 0; i < nPop; i++) {
		    rand_tmp = arma::randperm(nPop - 1, n_rand);
		    rand_tmp.elem(arma::find(rand_tmp >= i)) += 1;
		    rand_pop.col(i) = rand_tmp;
		  }
		  // the update formula of DE/rand/2
		  vMuta = pop.rows(rand_pop.row(0)) + sf*(pop.rows(rand_pop.row(1)) - pop.rows(rand_pop.row(2))) + sf*(pop.rows(rand_pop.row(3)) - pop.rows(rand_pop.row(4)));
		  // Generate new population
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
	  case 2: // DE/best/1
		{
		  arma::uword n_rand = 2;
		  arma::umat rand_pop(n_rand, nPop);
		  arma::uvec rand_tmp(n_rand);
		  for (int i = 0; i < nPop; i++) {
		    rand_tmp = arma::randperm(nPop - 1, n_rand);
		    rand_tmp.elem(arma::find(rand_tmp >= i)) += 1;
		    rand_pop.col(i) = rand_tmp;
		  }
		  // the update formula of DE/best/1
		  vMuta = GBmat + sf*(pop.rows(rand_pop.row(0)) - pop.rows(rand_pop.row(1)));
		  // Generate new population
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
	  case 3: // DE/best/2
		{
		  arma::uword n_rand = 4;
		  arma::umat rand_pop(n_rand, nPop);
		  arma::uvec rand_tmp(n_rand);
		  for (int i = 0; i < nPop; i++) {
		    rand_tmp = arma::randperm(nPop - 1, n_rand);
		    rand_tmp.elem(arma::find(rand_tmp >= i)) += 1;
		    rand_pop.col(i) = rand_tmp;
		  }
		  // the update formula of DE/best/2
		  arma::mat GBmat = repmat(GBest, nPop, 1);
		  vMuta = GBmat + sf*(pop.rows(rand_pop.row(0)) - pop.rows(rand_pop.row(1))) + sf*(pop.rows(rand_pop.row(2)) - pop.rows(rand_pop.row(3)));
		  // Generate new population
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
	  case 4: // DE/rand-to-best/1
		{
		  arma::uword n_rand = 2;
		  arma::umat rand_pop(n_rand, nPop);
		  arma::uvec rand_tmp(n_rand);
		  for (int i = 0; i < nPop; i++) {
		    rand_tmp = arma::randperm(nPop - 1, n_rand);
		    rand_tmp.elem(arma::find(rand_tmp >= i)) += 1;
		    rand_pop.col(i) = rand_tmp;
		  }
		  // the update formula of DE/rand-to-best/1
		  vMuta = pop + sf*(GBmat - pop) + sf*(pop.rows(rand_pop.row(0)) - pop.rows(rand_pop.row(1)));
		  // Generate new population
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
	  case 5: // DE/rand-to-best/2
		{
		  arma::uword n_rand = 4;
		  arma::umat rand_pop(n_rand, nPop);
		  arma::uvec rand_tmp(n_rand);
		  for (int i = 0; i < nPop; i++) {
		    rand_tmp = arma::randperm(nPop - 1, n_rand);
		    rand_tmp.elem(arma::find(rand_tmp >= i)) += 1;
		    rand_pop.col(i) = rand_tmp;
		  }
		  // the update formula of DE/rand-to-best/2
		  vMuta = pop + sf*(GBmat - pop) + sf*(pop.rows(rand_pop.row(0)) - pop.rows(rand_pop.row(1))) + sf*(pop.rows(rand_pop.row(2)) - pop.rows(rand_pop.row(3)));
		  // Generate new population
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
	PutRNGstate();
}
