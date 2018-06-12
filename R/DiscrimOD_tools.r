#' @export
designV2M <- function(v_design, D_INFO) {
	n <- D_INFO$nSupp
	d <- D_INFO$dSupp
	m_design <- switch (D_INFO$D_TYPE, 
		"exact" = {
			as.vector(t(v_design))
		},
		"approx" = {
			tmp <- as.vector(t(v_design[,1:d]))		
			wt <- sqrt(as.vector(v_design[,ncol(v_design)]))
			ang <- numeric(n-1)
			for (i in 1:(n-1)) {
				ang[i] <- acos(wt[i]/sqrt(sum((wt[i:n]^2))))
				if (i < (n-1)) {
			    ang[i] <- pi - ang[i]
			  }
			}
			c(tmp, ang)
		}
	)
	return(m_design)
}

#' @export
designM2V <- function(m_design, D_INFO) {
	n <- D_INFO$nSupp
	d <- D_INFO$dSupp
	v_design <- switch (D_INFO$D_TYPE, 
		"exact" = {
			tmp <- matrix(m_design, n, d, byrow = TRUE)
			tmp <- as.matrix(tmp[do.call(order,as.data.frame(round(tmp, 4))),])
			dimnames(tmp) <- list(paste0("obs_", 1:n), paste0("dim_", 1:d))
			tmp
		},
		"approx" = {
			# weight of support points
			ang <- m_design[(n*d + 1):length(m_design)]
			wcumsin <- wcos <- numeric(n)
		 	wcumsin[1] <- 1; wcumsin[2:n] <- cumprod(sin(ang))
		 	wcos[1:(n-1)] <- cos(ang); wcos[n] <- 1
			wt <-	(wcumsin*wcos)^2
			tmp <- cbind(matrix(m_design[1:(n*d)], n, d, byrow = TRUE), wt)
			tmp <- as.matrix(tmp[do.call(order, as.data.frame(round(tmp, 4))),])
			dimnames(tmp) <- list(paste0("obs_", 1:n), c(paste0("dim_", 1:d), "weight"))
			tmp
		},
		"maxmin_eqv_wt" = {
			n_model <- length(m_design) + 1
			ang <- m_design
			wcumsin <- wcos <- numeric(n_model)
			wcumsin[1] <- 1; wcumsin[2:n_model] <- cumprod(sin(ang))
		 	wcos[1:(n_model-1)] <- cos(ang); wcos[n_model] <- 1
			v_design <-	(wcumsin*wcos)^2
		}
	)
	return(v_design)
}

#' @export
algInfoUpdate <- function(D_INFO) {
	
	D_SWARM <- switch(D_INFO$D_TYPE, 
		"approx" = 
			list(UB = matrix(c(rep(D_INFO$dsUpper, D_INFO$nSupp), rep(pi  , D_INFO$nSupp - 2), pi/2), nrow = 1),
					 LB = matrix(c(rep(D_INFO$dsLower, D_INFO$nSupp), rep(pi/2, D_INFO$nSupp - 2),  0.0), nrow = 1)
					),
		"exact" = 
			list(UB = matrix(c(rep(D_INFO$dsUpper, D_INFO$nSupp)), nrow = 1),
				   LB = matrix(c(rep(D_INFO$dsLower, D_INFO$nSupp)), nrow = 1)
					),
		"multiexact" = 
			list(UB = matrix(c(rep(D_INFO$dsUpper, D_INFO$nSubj*D_INFO$nSupp), rep(pi  , D_INFO$nSubj - 2), pi/2), nrow = 1),
		 			 LB = matrix(c(rep(D_INFO$dsLower, D_INFO$nSubj*D_INFO$nSupp), rep(pi/2, D_INFO$nSubj - 2),  0.0), nrow = 1)
					),
		"2Dgrid" = 
			list(UB = matrix(c(rep(D_INFO$dsUpper[1], D_INFO$nGrid1), rep(D_INFO$dsUpper[2], D_INFO$nGrid2)), nrow = 1),
					 LB = matrix(c(rep(D_INFO$dsLower[1], D_INFO$nGrid1), rep(D_INFO$dsLower[2], D_INFO$nGrid2)), nrow = 1)
					),
		"maxmin_eqv_wt" = 
			list(UB = matrix(rep(pi/2, D_INFO$N_PAIR - 1), nrow = 1),
					 LB = matrix(rep(0.0,  D_INFO$N_PAIR - 1), nrow = 1)
					)			
	)
	D_SWARM
}
