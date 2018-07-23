################################################################
# Internal functions to be used in the BUSgibbs function
################################################################

# A function for sampling from the Dirichlet distribution
.Dirichlet <- function(alpha_vec){
    #input is a parameter vector
    #output is one sample from Dirichlet distribution
    
    n <- length(alpha_vec)
    gamma_rvs <- vapply(seq_len(n), function(i){
        rgamma(1, shape = alpha_vec[i], rate = 1) }, FUN.VALUE=1)
    return(gamma_rvs/sum(gamma_rvs))
}


# A function for sampling from the Inverse Gamma distribution
.Inv_Gamma <- function(a,b){
    #input is shape a and scale b of Inv-Gamma distribution
    #output is one sample from the inverse gamma distribution
    tmp <- rgamma(1,shape = a, rate = b)
    return(1/tmp)
}




################################################################
# Internal functions to be used in the BUSgibbs function
# for Gibbs sampler
################################################################


#Throughout the code, G is the gene number, K is the subtype
# number, B is the batch number, and n_vec is the sample size
# vector


# Full Conditional Functions

#Update the differentially expressed (DE) indicator proportion
.update_DE_prop <- function(L_t, a_p, b_p, G, K){
    #a_p and b_p are the hyperparameters for the DE proportion parameter
    s <- sum(L_t)
    tmp <- rbeta(1, shape1 = s + a_p, shape2 = G*(K-1) - s + b_p)
    return(tmp)
}

#Update the subtype proportions
.update_pi <- function(Z_t, alpha_par, B, K){
    #update the subtype proportion
    #Z_t is a list containing B n_vec[b]-component vector where the entry
    #Z_t[[b]][j] stands for the subtype indicator of jth subject in bth batch
    #alpha_par is the hyperparameter vector in the Dirichlet prior
    #output is a B by K matrix
    tmp <- matrix(NA, B, K)
    for(b in seq_len(B)){
        prob <- vapply(seq_len(K), function(k) sum(Z_t[[b]]==k), FUN.VALUE=1 )
        tmp[b, ] <- .Dirichlet(prob + alpha_par)
    }
    return(tmp)
}


#Update the DE indicators
.update_L <- function(mu_t, DE_prop_t, tau_mu_zero_t, tau_mu_one_t){
    #update L_{gk} for k >= 2, L_{gk} = 1 if the expression level of
    #subtype k's gene g is differentially expressed compared to that
    #of subtype 1's gene g.
    #output is a G by K-1 matrix
    args <- list("mu_t" = as.numeric(mu_t), "mu_t_dim" = as.integer(dim(mu_t)),
                 "DE_prop_t" = as.numeric(DE_prop_t),
                 "tau_mu_zero_t" = as.numeric(tau_mu_zero_t),
                 "tau_mu_one_t"=as.numeric(tau_mu_one_t))
    L <- .Call("update_L_c", args)
    return(L)
}

#Update the standard deviation for the spike part
.tau_mu_zero_update <- function(L_t, mu_t, a_tau0, b_tau0){
    ind <- (L_t == 0)
    s1 <- sum(ind)
    mu_t_tmp <- mu_t[ ,-1]
    s2 <- sum(mu_t_tmp[ind]^2)
    return(sqrt(.Inv_Gamma(a_tau0 + 1/2*s1, b_tau0 + 1/2*s2)))
}

#Update the subtypes for samples
.update_Z_v2 <- function( Z_t, alpha_t, mu_t,gamma_t,sigma_sq_t, pi_t,
                          Y, B, n_vec, K){
    #update Z by Metropolis-Hasting step
    #Z_t is a list containing B n_vec[b]-component vector where the entry
    #Z_t[[b]][j] stands for the subtype indicator of jth subject in bth batch
    #mu_t is the G by K subtype effect matrix
    #gamma_t is the B by G batch effect matrix
    #sigma_sq_t is the B by G batch effect (variance) matrix
    #pi_t is a B by K matrix
    #Y is the gene expression data list
    #output is a list containing B n_vec[b]-component vector
    
    Z_tmp <- Z_t
    for(b in seq_len(B)){
        for(j in seq_len(n_vec[b])){
            #get a proposal
            proposal <- sample(seq_len(K),1)
            
            #calculate posterior ratio in log scale
            aa <- sum( -(Y[[b]][j,] - alpha_t - mu_t[,proposal] -
                       gamma_t[b,])^2/(2*sigma_sq_t[b,]) )
                       
            bb <- sum( -(Y[[b]][j,] - alpha_t - mu_t[,Z_tmp[[b]][j]] -
                       gamma_t[b,])^2/(2*sigma_sq_t[b,]) )
            
            #calculate the acceptance prob
            prob <- min(exp(aa - bb)*pi_t[b, proposal]/pi_t[b, Z_tmp[[b]][j]] ,
                        1)
                        
            tmp <- runif(1)
            
            #accept or reject the proposal
            if(tmp <= prob){
                Z_tmp[[b]][j] <- proposal
            }
        }
    }
    return(Z_tmp)
}

#Update the baseline expression level
.update_alpha <- function(mu_t, gamma_t, Z_t, sigma_sq_t, tau_alpha,
                          eta_alpha, Y, n_vec){
    #update baseline expression level
    #mu_t is the G by K subtype effect matrix
    #gamma_t is the B by G batch effect matrix
    #Z_t is a list containing B n_vec[b]-component vector where the entry
    #Z_t[[b]][j] stands for the subtype indicator of jth subject in bth batch
    #sigma_sq_t is the B by G batch effect (variance) matrix
    #output is a G dimensional vector
    args <- list("Y" = Y, "n_vec" = as.integer(n_vec),
                 "mu_t" = as.numeric(mu_t),
                 "mu_t_dim" = as.integer(dim(mu_t)),
                 "gamma_t" = as.numeric(gamma_t),
                 "gamma_t_dim" = as.integer(dim(gamma_t)),
                 "Z_t" = Z_t, "sigma_sq_t" =  as.numeric(sigma_sq_t),
                 "tau_alpha" = as.numeric(tau_alpha),
                 "eta_alpha" = as.numeric(eta_alpha))
    
    alpha <- .Call("update_alpha_c", args)
    return(alpha)
}


#Update the subtype effects
.update_mu <- function(L_t, alpha_t, gamma_t, Z_t, sigma_sq_t, tau_mu_zero_t,
                       tau_mu_one_t, Y, B, n_vec, G, K){
    #update subtype effect
    #L_t is the G by K-1 matrix
    #alpha_t is the G dimensional vector
    #gamma_t is the B by G batch effect matrix
    #Z_t is a list containing B n_vec[b]-component vector where the entry
    #Z_t[[b]][j] stands for the subtype indicator of jth subject in bth batch
    #sigma_sq_t is the B by G batch effect (variance) matrix
    #output is a G by K matrix
    args <- list("Y" = Y, "L_t" = as.integer(L_t),
                 "alpha_t" = as.numeric(alpha_t),
                 "gamma_t" = as.numeric(gamma_t),
                 "Z_t" = Z_t, "sigma_sq_t" = as.numeric(sigma_sq_t),
                 "tau_mu_zero_t" = as.numeric(tau_mu_zero_t),
                 "tau_mu_one_t" = as.numeric(tau_mu_one_t),
                 "B" = as.integer(B), "n_vec" = as.integer(n_vec),
                 "G" = as.integer(G), "K" = as.integer(K))
    mu <- .Call("update_mu_c", args)
    return(mu)
}

#Update the location batch effects
.update_gamma <- function( alpha_t, mu_t, sigma_sq_t, Z_t, tau_gamma,
                           Y, B, n_vec, G){
    #update batch effect
    #alpha_t is the G dimensional vector
    #mu_t is a G by K matrix
    #sigma_sq_t is the B by G batch effect (variance) matrix
    #Z_t is a list containing B n_vec[b]-component vector where the entry
    #Z_t[[b]][j] stands for the subtype indicator of jth subject in bth batch
    #output is a B by G matrix
    
    args <- list("Y" = Y, "alpha_t" = as.numeric(alpha_t),
                 "mu_t" = as.numeric(mu_t),
                 "Z_t" = Z_t, "sigma_sq_t" = as.numeric(sigma_sq_t),
                 "B" = as.integer(B), "n_vec" = as.integer(n_vec),
                 "G" = as.integer(G),
                 "tau_gamma" = as.numeric(tau_gamma))
    
    gamma <- .Call("update_gamma_c", args)
    return(gamma)
}

#Update the variances within each batch
.update_sigma_sq <- function(alpha_t, mu_t, Z_t, gamma_t, a_inv_gamma,
                             b_inv_gamma, Y, B, n_vec, G){
    #update variance
    #alpha_t is the G dimensional vector
    #mu_t is the G by K subtype effect matrix
    #Z_t is a list containing B n_vec[b]-component vector where the entry
    #Z_t[[b]][j] stands for the subtype indicator of jth subject in bth batch
    #gamma_t is the B by G batch effect matrix
    #output is a B by G matrix
    
    args <- list("Y" = Y, "alpha_t" = as.numeric(alpha_t),
                 "mu_t" = as.numeric(mu_t),
                 "Z_t" = Z_t, "gamma_t" = as.numeric(gamma_t),
                 "B" = as.integer(B), "n_vec" = as.integer(n_vec),
                 "G" = as.integer(G), "a_inv_gamma" = as.numeric(a_inv_gamma),
                 "b_inv_gamma" = as.numeric(b_inv_gamma))
    
    sigma_sq <- .Call("update_sigma_sq_c", args)
    return(sigma_sq)
}


#Function for calculating estimated Observed Log Likelihood
.observed_log_likelihood <- function(pi_t_post, alpha_t_post, mu_t_post,
                                     gamma_t_post, sigma_sq_t_post, Y, n_vec){
    args <- list("Y" = Y, "n_vec" = as.integer(n_vec),
                  "pi_t_post"=as.numeric(pi_t_post),
                  "alpha_t_post" = as.numeric(alpha_t_post),
                  "mu_t_post" = as.numeric(mu_t_post),
                  "mu_t_dim" = as.integer(dim(mu_t_post)),
                  "gamma_t_post" = as.numeric(gamma_t_post),
                  "gamma_t_dim" = as.integer(dim(gamma_t_post)),
                  "sigma_sq_t_post" = as.numeric(sigma_sq_t_post))
    
    ret_value <- .Call("observed_log_likelihood_c", args)
    return(ret_value)
}


################################################################
# The Main Function
################################################################

BUSgibbs <- function(Data, n.subtypes, n.iterations = 500, 
				     n.records = floor(n.iterations / 2),  hyperparameters = 
	                 c(1, sqrt(5), sqrt(5), 2, 2, 1, 2, 0.005, 1, 3, 10),
                    showIteration = TRUE ){

	# Gibbs sampler
	Y <- list() #input data Y
	
	if(is(Data, "SummarizedExperiment")){
        #The input data format is SummarizedExperiment
		batch_info <- colData(Data)$Batch
		B <- length(unique(batch_info)) #batch number
		G <- nrow(assays(Data)$GE_matr)      #gene number
		n_vec <- NULL                   #sample size vector 
		for(b in seq_len(B)){
			n_vec <- c(n_vec, sum(batch_info == b))
			Y[[b]] <- t(assays(Data)$GE[ , batch_info == b]) #input data
		}
	}else if(is(Data, "list")){         #The input data format is list
		Y0 <- Data
		B <- length(Y0)                 #batch number    
		for(b in seq_len(B)){           #input data
			Y[[b]] <- t(Y0[[b]])
		}
		n_vec <- NULL                   #sample size vector 
		for(b in seq_len(B)){
			n_vec <- c(n_vec, dim(Y[[b]])[1] )
		}
		
		G_vec <- NULL                   #gene number vector
		for(b in seq_len(B)){
			G_vec <- c(G_vec, dim(Y[[b]])[2] )
		}

		if(length(unique(G_vec)) == 1){ #consistent gene numbers
			G <- G_vec[1]
		}else stop("The gene numbers across batches must be the same.\n")

		if(sum(G_vec <= n_vec) > 0){
            #the gene number must be greater than
            #the sample size in each batch
			stop("The gene number must be great than the sample size.\n")
		}
	}else{
		stop(paste0("Data must be either a \"SummarizedExperiment\" object",
                    " or a \"list\" object!\n"))
	}
	
	if(B < 2){
		stop("The batch number must be greater than one.\n")
	}


	K <- n.subtypes

	if(sum(K > n_vec) > 0){
		stop(paste0("The sample size in any batch must be greater",
                    " than the assumed subtype number.\n"))
	}

	###specify hyperparameters

	eta_alpha <- hyperparameters[1]
	tau_alpha <- hyperparameters[2]
	tau_gamma <- hyperparameters[3]
	alpha_par <- hyperparameters[4]
	a_inv_gamma <- hyperparameters[5]
	b_inv_gamma <- hyperparameters[6]
	a_tau0 <- hyperparameters[7]
	b_tau0 <- hyperparameters[8]
	a_p <- hyperparameters[9]
	b_p <- hyperparameters[10]

	###single chain

	###set initial values
	DE_prop_t <- runif(1, 0, 1/2)
	tau_mu_zero_t <- 1
	tau_mu_one_t <- hyperparameters[11]  #tau_mu_one_t are always fixed
	pi_t <- array(NA, dim = c( B, K))

	###initialize Z_t
	Z_t <- list()

	for(b in seq_len(B)){
		repeat{
			Z_t[[b]] <- sample(seq_len(K), n_vec[b], replace = TRUE) 
			if(length(unique(Z_t[[b]])) == K){
				break
			}
		}
	}


	###set data-dependent initial values 
	###to alpha_t

	raw_Means<- array(NA, dim = c(B, K, G))

    #Note: the product of the batch number B and the subtype number K
    #is usually less than 100. The following nested for loop cost
    #little time.
    
	for(b in seq_len(B)){
		for(k in seq_len(K)){
            sumzt_bk <- sum(Z_t[[b]] == k)
			if(sumzt_bk > 1){
				raw_Means[b,k,] <- colMeans(Y[[b]][Z_t[[b]]==k,])
            }else{
                raw_Means[b,k,] <- Y[[b]][Z_t[[b]]==k,]
            }
		}
	}

	alpha_t <- raw_Means[1,1, ]


	###to mu_t
	mu_t <- array(NA, dim = c( G, K))
	mu_t[ ,1] <- 0


	for(k in 2:K){
		mu_t[ ,k] <- raw_Means[1,k, ] - alpha_t
	}


	###to L_t
	L_t <- .update_L(mu_t, DE_prop_t, tau_mu_zero_t, tau_mu_one_t)

	###to tau_mu_zero_t
	tau_mu_zero_t <- .tau_mu_zero_update(L_t, mu_t, a_tau0, b_tau0)


	###to gamma_t
	gamma_t <- array(0,dim = c(B, G))

	for(b in 2:B){
		gamma_t[b, ] <- colMeans(raw_Means[b,,]-raw_Means[1,,]) 
	}



	###to sigma_sq_t

	sigma_sq_t <- .update_sigma_sq(alpha_t, mu_t, Z_t, gamma_t,
				a_inv_gamma, b_inv_gamma, Y, B, n_vec, G)


	# record samples after burn-in period

	T_iter <- n.iterations #total number of iterations
	num_rec <- n.records   #how many last iterations to be used
                           #to infer parameters


	### variables used to record samples
	Z_t_record <- list()
	for(b in seq_len(B)){
		Z_t_record[[b]] <- matrix(NA, n_vec[b], num_rec)
	}


	DE_prop_t_record <- rep(NA, num_rec)
	alpha_t_record <- array(NA, dim = c(G, num_rec))
	L_t_record <- array(NA, dim = c(G, K-1, num_rec))
	pi_t_record <- array(NA, dim = c( B, K, num_rec))
	mu_t_record <- array(NA, dim = c( G, K, num_rec))
	gamma_t_record <- array(NA, dim = c( B, G, num_rec))
	sigma_sq_t_record <- array(NA, dim = c( B, G, num_rec))
	tau_mu_zero_t_record <- rep(NA, num_rec)

	t1 <- Sys.time()
    
    #Gibbs sampler begins here
	message("  running the Gibbs sampler ...\n")
	for(t in seq_len(T_iter)){
        #sequetially update each parameter
		pi_t <- .update_pi(Z_t, alpha_par, B, K)
		DE_prop_t <- .update_DE_prop(L_t, a_p, b_p, G, K)
		Z_t <- .update_Z_v2(Z_t, alpha_t, mu_t, gamma_t,sigma_sq_t,pi_t,
                            Y, B, n_vec, K)
		L_t <- .update_L(mu_t, DE_prop_t, tau_mu_zero_t, tau_mu_one_t)
		tau_mu_zero_t <- .tau_mu_zero_update(L_t, mu_t, a_tau0, b_tau0)
		alpha_t <- .update_alpha(mu_t, gamma_t, Z_t, sigma_sq_t, tau_alpha,
                                 eta_alpha, Y, n_vec)
		mu_t <- .update_mu(L_t, alpha_t, gamma_t, Z_t, sigma_sq_t,
                           tau_mu_zero_t,  tau_mu_one_t, Y, B, n_vec, G, K)
		gamma_t <- .update_gamma(alpha_t, mu_t, sigma_sq_t, Z_t, tau_gamma,
                                 Y, B, n_vec, G)
		sigma_sq_t <- .update_sigma_sq(alpha_t, mu_t, Z_t, gamma_t,
                                  a_inv_gamma, b_inv_gamma, Y, B, n_vec, G)
	
		if(showIteration == TRUE){
			message(c("  Iteration ", t, "\n") )
		}

        #when the iteration number is large enough,
        #the samples are stored
		if(t > T_iter - num_rec ){
			DE_prop_t_record[t - (T_iter - num_rec)] <- DE_prop_t
			tau_mu_zero_t_record[t - (T_iter - num_rec)] <- tau_mu_zero_t
			alpha_t_record[ , t - (T_iter - num_rec)] <- alpha_t
			pi_t_record[ , ,t - (T_iter - num_rec)] <- pi_t
			L_t_record[ , ,t - (T_iter - num_rec)] <- L_t
			mu_t_record[ , ,t - (T_iter - num_rec)] <- mu_t
			gamma_t_record[ , ,t - (T_iter - num_rec)] <- gamma_t
			sigma_sq_t_record[ , ,t - (T_iter - num_rec)] <- sigma_sq_t
			for(b in seq_len(B)){
				Z_t_record[[b]][ ,t - (T_iter - num_rec)] <- Z_t[[b]]
			}
		
		}

	} 

	t2 <- Sys.time()

	output <- list()

	message(paste0("  The Gibbs sampler takes: ", round(difftime( t2, t1,
                units = "mins"), 3), " mins", "\n"))
	###Posterior samples of DE indicators L

	message("  calculating posterior means and posterior modes...\n")

	output[[1]] <- L_t_record

	###Posterior mode of Z
	Z_post <- list()

	for(b in seq_len(B)){
		Z_post[[b]] <- vapply(seq_len(n_vec[b]), function(j){

							temp <- table(Z_t_record[[b]][j,])
							as.numeric(names(temp)[which.max(temp)])
						}, FUN.VALUE=1)

	}

	output[[2]] <- Z_post

	###Posterior mean of parameters

	pi_t_post <- matrix(NA, B, K)
	alpha_t_post <- rep(NA,G)
	tau_mu_zero_t_post <- NA
	DE_prop_t_post <- NA
	mu_t_post <- matrix(NA, G, K)
	gamma_t_post <- matrix(NA, B, G)
	sigma_sq_t_post <- matrix(NA, B, G)

	tau_mu_zero_t_post <- mean(tau_mu_zero_t_record)
	DE_prop_t_post <- mean(DE_prop_t_record)


	for(b in seq_len(B)){
		pi_t_post[b, ] <- rowMeans(pi_t_record[b,seq_len(K), ])
		gamma_t_post[b, ] <- rowMeans(gamma_t_record[b,seq_len(G), ])
		sigma_sq_t_post[b, ] <- rowMeans(sigma_sq_t_record[b,seq_len(G), ])
	}


	for(g in seq_len(G)){
		alpha_t_post[g] <- mean(alpha_t_record[g, ])
		mu_t_post[g, ] <- rowMeans(mu_t_record[g,seq_len(K), ])
	}


	output[[3]] <- tau_mu_zero_t_post 

	output[[4]] <- DE_prop_t_post 

	output[[5]] <- pi_t_post

	output[[6]] <- alpha_t_post

	#transpose begin
	output[[7]] <- aperm(gamma_t_record, c(2,1,3))
	output[[8]] <- t(gamma_t_post) 

	output[[9]] <- aperm(sigma_sq_t_record, c(2,1,3))
	output[[10]] <- t(sigma_sq_t_post)
	#transpose end

	output[[11]] <- mu_t_record
	output[[12]] <- mu_t_post

	message("  calculating BIC...\n")
	BIC <- (-2)*.observed_log_likelihood(pi_t_post, alpha_t_post, mu_t_post,
                                gamma_t_post, sigma_sq_t_post, Y, n_vec) +
                 ((2*B-1)*G + G*K)*log(sum(n_vec)*G)

	output[[13]] <- BIC


	names(output) <- c("L_PosterSamp", "Subtypes", "tau_mu_zero", "p", "pi",
                       "alpha", "gamma_PosterSamp", "gamma",
			           "sigma_sq_PosterSamp", "sigma_sq", "mu_PosterSamp",
                       "mu", "BIC")

	class(output) <- "BUSfits"
	return(output)

}

############################################################################
#Estimate intrinsic gene indicators
############################################################################
	
#calculate posterior probability of being differentially expressed 
#for gene g in subtype k (k>=2) compared to subtype 1	
postprob_DE <- function(BUSfits){
	if(!is(BUSfits, "BUSfits")){
		stop("BUSfits should be in the \"BUSfits\" class.\n")
	}
	L_PosterSamp <- BUSfits$L_PosterSamp
	G <- dim(L_PosterSamp)[1]
	K <- dim(L_PosterSamp)[2] + 1
	num_rec <- dim(L_PosterSamp)[3]
	PPI <- NULL
	for(k in seq_len(K-1)){
		temp <- vapply(seq_len(G), function(j){
					sum(L_PosterSamp[j,k, ] == 1) / num_rec
				}, FUN.VALUE=1)
		PPI <- cbind(PPI, temp)
	}
	message("Showing the posterior probability of being differentially",
            "expressed\n")
	message(paste0("               for gene g in subtype k (k>=2) compared",
               " to subtype 1.\n\n"))
	message("The output format is a matrix.\n")
	message(paste0("Each row represents a gene, and each column corresponds",
               " to a subtype.\n"))
	return(PPI)
}



#Internal functions used in the next function
#"postprob_DE_thr_fun"
.postprob_DE2 <- function(BUSfits){
    if(!is(BUSfits, "BUSfits")){
        stop("BUSfits should be in the \"BUSfits\" class.\n")
    }
    L_PosterSamp <- BUSfits$L_PosterSamp
    G <- dim(L_PosterSamp)[1]
    K <- dim(L_PosterSamp)[2] + 1
    num_rec <- dim(L_PosterSamp)[3]
    PPI <- NULL
    for(k in seq_len(K-1)){
        temp <- vapply(seq_len(G), function(j){
            sum(L_PosterSamp[j,k, ] == 1) / num_rec
        }, FUN.VALUE=1)
        PPI <- cbind(PPI, temp)
    }
    
    return(PPI)
}

.fdrDEindicator <- function(L_PosterSamp, kappa){
    G <- dim(L_PosterSamp)[1]
    K <- dim(L_PosterSamp)[2] + 1
    num_rec <- dim(L_PosterSamp)[3]
    
    args <- list("G"=as.integer(G), "K"=as.integer(K),
                 "num_rec"=as.integer(num_rec), "kappa"=as.numeric(kappa),
                 "L_PosterSamp"=as.integer(L_PosterSamp))
    fdr <- .Call("fdrDEindicator_c", args)
    
    return(fdr)
}

#calculate the DE posterior probability threshold
postprob_DE_thr_fun <- function(BUSfits, fdr_threshold=0.1){
	#find posterior probability threshold to control FDR
	if(!is(BUSfits, "BUSfits")){
		stop("BUSfits should be in the \"BUSfits\" class.\n")
	}
	L_PosterSamp <- BUSfits$L_PosterSamp
	#L_PosterSamp is the posterior samples of the intrinsic gene indicators
	#alpha (default is 0.1) is the threshould of fdr we want to control	

	kappa_fdr_matr <- NULL
	kappa_set <- rev(1 - sort(unique(c(.postprob_DE2(BUSfits)))))
	for(kappa in kappa_set){
		fdr <- .fdrDEindicator(L_PosterSamp, kappa=kappa)
		kappa_fdr_matr <- rbind(kappa_fdr_matr, c(kappa, fdr))
		if(fdr > fdr_threshold){
			break
		}
	}
	ind <- which(kappa_fdr_matr[ ,2] <= fdr_threshold)
	ind2 <- which.max(kappa_fdr_matr[ind,2])
	ind3 <- ind[ind2] # the index that has the maximum fdr
                      # but less than fdr_threshold
	message(c("Posterior probability threshold = ",
          1-as.numeric(kappa_fdr_matr[ind3,1]),"\n"))
	message("The output is a scalar.\n")
	return(1-as.numeric(kappa_fdr_matr[ind3,1])) #1-kappa
}

#Estimate intrinsic gene indicators
estimate_IG_indicators <- function(BUSfits, postprob_DE_threshold = 0.5){
	if(!is(BUSfits, "BUSfits")){
		stop("BUSfits should be in the \"BUSfits\" class.\n")
	}
	L_PosterSamp <- BUSfits$L_PosterSamp
	G <- dim(L_PosterSamp)[1]
	K <- dim(L_PosterSamp)[2] + 1
	num_rec <- dim(L_PosterSamp)[3]
	PPI <- NULL
	for(k in seq_len(K-1)){
		temp <- vapply(seq_len(G), function(j){
					sum(L_PosterSamp[j,k, ] == 1) / num_rec
				}, FUN.VALUE=1)
		PPI <- cbind(PPI, temp)
	}
	EstL <- PPI
	EstL[PPI >= postprob_DE_threshold] <- 1
	EstL[PPI < postprob_DE_threshold] <- 0
	message("The output format is a matrix.\n")
	message(paste0("Each row represents a gene, and each column",
               " corresponds to a subtype from 2 to K\n"))
	return(EstL)
}

#Intrinsic gene index
IG_index <- function(EstIGindicators){
	ind <- which(rowSums(EstIGindicators) > 0)
	message(c(length(ind), "intrinsic genes are found.\n"))
	message("The output format is a vector showing the intrinsic gene",
            " indices.\n")
	return(ind)
}

################################################################
# Adjusted Values
################################################################
#calculate the adjusted gene expression values
adjusted_values <- function(BUSfits, original_data){
	if(!is(BUSfits, "BUSfits")){
		stop("BUSfits should be in the \"BUSfits\" class.\n")
	}
	Y0 <- original_data
	B <- length(Y0)
	#transpose
	Y <- list()
	for(b in seq_len(B)){
		Y[[b]] <- t(Y0[[b]])
	}
	n_vec <- NULL #sample size 
	for(b in seq_len(B)){
		n_vec <- c(n_vec, nrow(Y[[b]]))
	}
	G <- ncol(Y[[1]]) #gene number

	Y_correct <- list()
	Y_correct[[1]] <- t(Y[[1]])
	for(b in 2:B){
		matr <- matrix(NA, n_vec[b],G)
		for(j in seq_len(n_vec[b])){
			matr[j, ] <- BUSfits$alpha +
                         BUSfits$mu[, BUSfits$Subtypes[[b]][j]] +
						 (Y[[b]][j,] - BUSfits$alpha -
                             BUSfits$mu[,BUSfits$Subtypes[[b]][j]] -
                             BUSfits$gamma[,b]) /
                          sqrt(BUSfits$sigma_sq[,b]/BUSfits$sigma_sq[,1])
		}
		Y_correct[[b]] <- t(matr)
	} 
	message("The output format is a list with length equal to",
            " the batch number.\n")
	message("Each element of the list is the adjusted gene",
            " expression matrix.\n")
	message(paste0("In the matrix, each row represents a gene,",
               " and each column corresponds to a sample.\n"))
	return(Y_correct)
}

##############################################################################
# Useful Outputs from BUSfits
##############################################################################
#obtain the subtype indicators for samples
Subtypes <- function(BUSfits){
	if(!is(BUSfits, "BUSfits")){
		stop("BUSfits should be in the \"BUSfits\" class.\n")
	}
	B <- length(BUSfits$Subtypes)
	for(b in seq_len(B)){
		message(c("Batch ", b, " samples' subtype indicators: ",
            BUSfits$Subtypes[[b]][1],
			BUSfits$Subtypes[[b]][2], BUSfits$Subtypes[[b]][3], "... ...\n"))
	}
	message("The output format is a list with length",
            " equal to the batch number.\n")
	message(paste0("Each element of the list is a subtype indicator vector in",
               " that batch.\n"))
	return(BUSfits$Subtypes)
}

#obtain the baseline expression values
baseline_expression_values <- function(BUSfits){
	if(!is(BUSfits, "BUSfits")){
		stop("BUSfits should be in the \"BUSfits\" class.\n")
	}	
	message("The output format is a vector.\n")
	return(BUSfits$alpha)
}

#obtain the subtype effects
subtype_effects <- function(BUSfits){
	if(!is(BUSfits, "BUSfits")){
		stop("BUSfits should be in the \"BUSfits\" class.\n")
	}
	message("The output format is a matrix.\n")
	message(paste0("Each row represents a gene, and each column corresponds",
               " to a subtype.\n"))
	
	return(BUSfits$mu)
}

#obtain the location batch effects
location_batch_effects <- function(BUSfits){
	if(!is(BUSfits, "BUSfits")){
		stop("BUSfits should be in the \"BUSfits\" class.\n")
	}
	message("The output format is a matrix.\n")
	message("Each row represents a gene, and each column",
            " corresponds to a batch.\n")
	return(BUSfits$gamma)
}

#obtain the scale batch effects
scale_batch_effects <- function(BUSfits){
	if(!is(BUSfits, "BUSfits")){
		stop("BUSfits should be in the \"BUSfits\" class.\n")
	}
	G <- nrow(BUSfits$sigma_sq)
	B <- ncol(BUSfits$sigma_sq)
	message("The output format is a matrix.\n")
	message(paste0("Each row represents a gene, and each column corresponds",
               " to a batch.\n"))
	tmp <- sqrt(BUSfits$sigma_sq / BUSfits$sigma_sq[,1])
	return(tmp)
}

#obtain the BIC score
BIC_BUS <- function(BUSfits){
	message("BIC is ", BUSfits$BIC, "\n")
	message("The output is a scalar.\n")
	return(BUSfits$BIC)
}

##############################################################################
# print and summary
##############################################################################
#print BUSfits
print.BUSfits <- function(x, ...){
	BUSfits <- x
	G <- nrow(BUSfits$sigma_sq)
	B <- ncol(BUSfits$sigma_sq)
	cat("Subtype indicators:\n")
	for(b in seq_len(B)){
		cat(c("   Batch ", b, " samples' subtype indicators: ",
            BUSfits$Subtypes[[b]][1],
			BUSfits$Subtypes[[b]][2], BUSfits$Subtypes[[b]][3], "... ...\n"))
	}
	cat("\n")
	cat("Estimated location batch effects:\n")
	for(b in seq_len(B)){
		cat(c("   Batch ", b, " location batch effects are: ",
            BUSfits$gamma[1,b],
			BUSfits$gamma[2,b], BUSfits$gamma[3,b], "... ...\n"))
	}
	cat("\n")
	cat("Estimated scale batch effects:\n")
	for(b in seq_len(B)){
		cat(c("   Batch ", b, " scale batch effects are: ",
            sqrt(BUSfits$sigma_sq[1,b]/BUSfits$sigma_sq[1,1]),
			sqrt(BUSfits$sigma_sq[2,b]/BUSfits$sigma_sq[2,1]), 
			sqrt(BUSfits$sigma_sq[3,b]/BUSfits$sigma_sq[3,1]), "... ...\n"))
	}
	cat("\n")
}

#summarize BUSfits
summary.BUSfits <- function(object, ...){
	BUSfits <- object
	G <- nrow(BUSfits$sigma_sq)
	B <- ncol(BUSfits$sigma_sq)
	K <- ncol(BUSfits$mu)
	num_records <- dim(BUSfits$L_PosterSamp)[3]
	cat(c("B = ", B, " batches\n"))
	cat(c("G = ", G, " genes\n"))
	cat(c("K = ", K, " subtypes\n"))
	cat(c("n.records = ", num_records," iterations are recorded.\n\n"))
	cat("BUSfits is an R list that contains the following main elements:\n\n")
	cat(paste0("   BUSfits$Subtypes : estimated subtype indicators,",
               " an R list with length B.\n"))
	cat(paste0("   BUSfits$pi : estimated subtype proportions across batches,",
               " a B by K matrix.\n"))
	cat(paste0("   BUSfits$alpha : estimated baseline expression levels,",
               " a vector with length G.\n"))
    cat(paste0("   BUSfits$gamma : estimated location batch",
               " effects a G by B matrix.\n"))
	cat("   BUSfits$mu : estimated subtype effects, a G by K matrix.\n")
	cat(paste0("   BUSfits$sigma_sq : estimated variances across batches,",
               " a G by B matrix.\n"))
	cat("   BUSfits$BIC : estimated BIC, a scalar.\n")
	cat(paste0("   BUSfits$L_PosterSamp : the posterior samples of the",
               " intrinsic gene indicators,\n"))
	cat("                          a G by K-1 by n.records array.\n")
	cat("   For more output values, please use \"?BUSgibbs\"\n")
	cat("\n")
}

###########################################################################
# Calculate EPSR factors
###########################################################################

#calculate EPSR factors for subtype effects
calculate_EPSR_mu <- function(mu_PosterSamp_chain1, mu_PosterSamp_chain2){
		message("  calculating EPSR factors ...\n")
		G <- dim(mu_PosterSamp_chain1)[1]
		K <- dim(mu_PosterSamp_chain1)[2]
		num_records <- dim(mu_PosterSamp_chain1)[3]
		num_chains <- 2
		if(num_records %% 2 != 0){
			num_records <- num_records - 1 
		}
		mu_t_record <- array(NA, dim = c(num_chains, G, K, num_records))
		mu_t_record[1,,,] <- mu_PosterSamp_chain1
		mu_t_record[2,,,] <- mu_PosterSamp_chain2
		
		
		mu_t_collect <- mu_t_record
        #within-sequence variable
		W_var_temp <- array(NA, dim = c(2*num_chains, G, K))
		W_var <- array(NA, dim = c(G, K))
        
        #between-sequence variable
		B_var <- array(NA, dim = c(G, K))
		mean_chains <- array(NA, dim = c(2*num_chains, G, K))
		EPSR <- array(NA, dim = c(G, K))

	
		for(k in 2:K){
			for(i in seq_len(num_chains)){
				temp <- vapply(seq_len(G), function(g){
						 c(  var(mu_t_collect[i, g, k,
                                              seq_len(num_records / 2)]),
							 var(mu_t_collect[i, g, k,
                                    (num_records / 2 + 1):(num_records)]),
							 mean(mu_t_collect[i, g, k,
                                     seq_len(num_records / 2)]),
							 mean(mu_t_collect[i, g, k,
                                    (num_records / 2 + 1):(num_records)]) )
						}, FUN.VALUE=rep(-1, 4))
				W_var_temp[2*(i-1)+1, ,k] <- temp[1, ]
				W_var_temp[2*(i-1)+2, ,k] <- temp[2, ]
				mean_chains[2*(i-1)+1, ,k] <- temp[3, ]
				mean_chains[2*(i-1)+2, ,k] <- temp[4, ]
			}
			B_var[ , k] <- num_records / 2 * apply(mean_chains[ , ,k] , 2, var)
			W_var[ , k] <- colMeans(W_var_temp[ , ,k])
			EPSR[ , k] <- sqrt( (num_records / 2 - 1 + B_var[ ,k] /
                                W_var[ ,k]) / (num_records / 2) )
		}
		return(EPSR)
	
}

#calculate EPSR factors for location batch effects
calculate_EPSR_gamma <- function(gamma_PosterSamp_chain1,
                                 gamma_PosterSamp_chain2){
		message("  calculating EPSR factors ...\n")
		B <- dim(gamma_PosterSamp_chain1)[1]
		G <- dim(gamma_PosterSamp_chain1)[2]
		num_records <- dim(gamma_PosterSamp_chain1)[3]
		num_chains <- 2
		if(num_records %% 2 != 0){
			num_records <- num_records - 1 
		}
		gamma_t_record <- array(NA, dim = c(num_chains, B, G, num_records))
		gamma_t_record[1,,,] <- gamma_PosterSamp_chain1
		gamma_t_record[2,,,] <- gamma_PosterSamp_chain2

		gamma_t_collect <- gamma_t_record
        #within-sequence variable
		W_var_temp <- array(NA, dim = c(2*num_chains, B, G))
		W_var <- array(NA, dim = c(B, G))
        #between-sequence variable
		B_var <- array(NA, dim = c(B, G))
		mean_chains <- array(NA, dim = c(2*num_chains,B, G))
		EPSR <- array(NA, dim = c(B, G))

		for(b in 2:B){
			for(i in seq_len(num_chains)){
				 temp <- vapply(seq_len(G), function(g){
						c( var(gamma_t_collect[i, b, g,
                                        seq_len(num_records / 2)]),
						   var(gamma_t_collect[i, b, g,
                                   (num_records / 2 + 1):(num_records)]),
				 		   mean(gamma_t_collect[i, b, g,
                                         seq_len(num_records / 2)]),
						   mean(gamma_t_collect[i, b, g,
                                    (num_records / 2 + 1):(num_records)]) )
						}, FUN.VALUE=rep(-1, 4))
				W_var_temp[2*(i-1)+1,b, ] <- temp[1,]
				W_var_temp[2*(i-1)+2,b, ] <- temp[2,]
				mean_chains[2*(i-1)+1,b, ] <- temp[3,]
				mean_chains[2*(i-1)+2,b, ] <- temp[4,]
			}
			B_var[b, ] <- num_records / 2 * apply(mean_chains[ ,b, ], 2, var)
			W_var[b, ] <- colMeans(W_var_temp[ ,b, ])
			EPSR[b, ] <- sqrt( (num_records / 2 - 1 + B_var[b, ] / W_var[b, ]) /
                         (num_records / 2) )
		
		}
		return(t(EPSR))
	
}

#calculate EPSR factors for variances in each batch
calculate_EPSR_sigma_sq <- function(sigma_sq_PosterSamp_chain1,
                                    sigma_sq_PosterSamp_chain2){
		message("  calculating EPSR factors ...\n")
		B <- dim(sigma_sq_PosterSamp_chain1)[1]
		G <- dim(sigma_sq_PosterSamp_chain1)[2]
		num_records <- dim(sigma_sq_PosterSamp_chain1)[3]
		num_chains <- 2
		
		if(num_records %% 2 != 0){
			num_records <- num_records - 1 
		}

		sigma_sq_t_record <- array(NA, dim = c(num_chains, B, G, num_records))
		sigma_sq_t_record[1,,,] <- sigma_sq_PosterSamp_chain1
		sigma_sq_t_record[2,,,] <- sigma_sq_PosterSamp_chain2
	
		sigma_sq_t_collect <- sigma_sq_t_record
        #within-sequence variable
		W_var_temp <- array(NA, dim = c(2*num_chains, B, G))
		W_var <- array(NA, dim = c(B, G))
        #between-sequence variable
		B_var <- array(NA, dim = c(B, G))
		mean_chains <- array(NA, dim = c(2*num_chains,B, G))
		EPSR <- array(NA, dim = c(B, G))

		for(b in seq_len(B)){
			for(i in seq_len(num_chains)){
				 temp <- vapply(seq_len(G), function(g){
						c( var(sigma_sq_t_collect[i, b, g,
                                 seq_len(num_records / 2)]),
						   var(sigma_sq_t_collect[i, b, g,
                                 (num_records / 2 + 1):(num_records)]),
				 		   mean(sigma_sq_t_collect[i, b, g,
                                  seq_len(num_records / 2)]),
						   mean(sigma_sq_t_collect[i, b, g,
                                 (num_records / 2 + 1):(num_records)]) )
						}, FUN.VALUE=rep(-1, 4))
				W_var_temp[2*(i-1)+1,b, ] <- temp[1,]
				W_var_temp[2*(i-1)+2,b, ] <- temp[2,]
				mean_chains[2*(i-1)+1,b, ] <- temp[3,]
				mean_chains[2*(i-1)+2,b, ] <- temp[4,]
			}
			B_var[b, ] <- num_records / 2 * apply(mean_chains[ ,b, ], 2, var)
			W_var[b, ] <- colMeans(W_var_temp[ ,b, ])
			EPSR[b, ] <- sqrt( (num_records / 2 - 1 + B_var[b, ] / W_var[b, ]) /
                               (num_records / 2) )
		
		}
		return(t(EPSR))
	
}



########################################################################
# Visualization
########################################################################
#visualize the gene expression data by stacking all gene expression matrices 
visualize_data <- function(Data, title_name="Heatmap", gene_ind_set,
                           color_key_range=seq(-0.5,8.5,1)){
	B <- length(Data)
	G <- nrow(Data[[1]])
	n_vec <- NULL
	for(b in seq_len(B)){
		n_vec <- c(n_vec, ncol(Data[[b]]))
	}
    #cell colors
	colfunc <- colorRampPalette(c("grey", "black"))
    #batch colors
	color_batch_func <- colorRampPalette(c("skyblue", "slateblue4"))

	color_batch <- color_batch_func(B)

	color_batch2 <- NULL

	for(b in seq_len(B)) {
		color_batch2 <- c(color_batch2, rep(color_batch[b], n_vec[b]))
	}
	Y <- NULL
	for(b in seq_len(B)){
		Y <- cbind(Y, Data[[b]])
	}
	Y1 <- Y[gene_ind_set, ]
	heatmap.2(Y1, col = colfunc(length(color_key_range)-1), scale = "none",
              key = TRUE, Colv=FALSE,Rowv=FALSE,
	          density.info = "none", trace = "none", dendrogram="none",
	          ylab = paste0(length(gene_ind_set), " Genes"),
              xlab = paste0(sum(n_vec)," Samples"),
	          main = title_name,labRow=FALSE,labCol=FALSE,
              breaks=color_key_range, symkey = FALSE,
              ColSideColors = color_batch2)
}

