\name{BUScorrect-package}
\alias{BUScorrect-package}
\alias{BUScorrect}
\docType{package}
\title{
\packageTitle{BUScorrect}
}
\description{
\packageDescription{BUScorrect}
}
\details{
BUS is capable of (a) correcting batch effects explicitly, (b) grouping samples that share similar characteristics into subtypes, (c) identifying features that distinguish subtypes, (d) allowing the number of subtypes to vary from batch to batch, (e) integrating batches from different platforms, and (f) enjoying a linear-order computation complexity. 
}
\author{
\packageAuthor{BUScorrect}

Maintainer: \packageMaintainer{BUScorrect}
}
\references{
Xiangyu Luo, Yingying Wei. Batch Effects Correction with Unknown Subtypes. Journal of the American Statistical Association. Accepted.
}

\examples{
###############################################################################
#Generate Simulation Data
###############################################################################
rm(list = ls(all = TRUE))  
set.seed(123)

B <- 3					
#total number of batches

K <- 3					
#total number of subtypes

G <- 3000					
#total number of genes

pi <- matrix(NA, B, K)			
# pi[b,k] stands for the proportion of the kth subtype in the bth batch

pi[1, ] <- c(0.2, 0.3, 0.5)
pi[2, ] <- c(0.4, 0.2, 0.4)
pi[3, ] <- c(0.3, 0.4, 0.3)

	
#total number of samples in each bacth.
n_vec <- rep(NA, B) 		
	
#n_vec[b] represents the total number of samples in batch b.
n_vec <- c(100, 110, 120)

#Data list
example_Data <- list()

#baseline expression level
alpha <- rep(2, G)



#subtype effect
mu <- matrix(NA, G, K)			
#subtype effect, mu[g,k] stands for the kth-subtype effect of gene g 

mu[ ,1] <- 0		
#the first subtype is taken as the baseline subtype		
#the subtype effect of subtype 1 is set to zero

mu[ ,2] <- c(rep(2,G/20), rep(0,G/20),rep(0, G-G/20-G/20))
mu[ ,3] <- c(rep(0,G/20), rep(2,G/20),rep(0, G-G/20-G/20)) 

#batch effect
gamma <- matrix(NA, B, G)		
#'location' batch effect of gene g in batch b 
gamma[1, ] <- 0		
#the first batch is taken as the reference batch without batch effects		
#the batch effect of batch 1 is set to zero
gamma[2, ] <- c(rep(3,G/5),rep(2,G/5),rep(1,G/5),
										rep(2,G/5),rep(3,G/5))
gamma[3, ] <- c(rep(1,G/5),rep(2,G/5),rep(3,G/5),
										rep(2,G/5),rep(1,G/5))

sigma_square <- matrix(NA, B,G)	
#sigma_square[b,g] denotes the error variance of gene g in batch b.

sigma_square[1,] <- rep(0.1, G)
sigma_square[2,] <- rep(0.2, G)
sigma_square[3,] <- rep(0.15, G)


Z <- list()			
#subtype indicator. Z[b,j] represents the subtype of sample j in batch b

Z[[1]] <- as.integer(c(rep(1,floor(pi[1,1]*n_vec[1])),rep(2,floor(pi[1,2]*n_vec[1])),
 rep(3,floor(pi[1,3]*n_vec[1]))))	
 
Z[[2]] <- as.integer(c(rep(1,floor(pi[2,1]*n_vec[2])),rep(2,floor(pi[2,2]*n_vec[2])), 
rep(3,floor(pi[2,3]*n_vec[2]))))	

Z[[3]] <- as.integer(c(rep(1,floor(pi[3,1]*n_vec[3])),rep(2,floor(pi[3,2]*n_vec[3])), 
rep(3,floor(pi[3,3]*n_vec[3]))))


for(b in 1:B){				#generate data 
	num <- n_vec[b]
	example_Data[[b]] <-  sapply(1:num, function(j){
						tmp <- alpha + mu[ ,Z[[b]][j]] + gamma[b, ] + 
							rnorm(G, sd = sqrt(sigma_square[b, ])) 
									
						tmp
					}) 

}

###############################################################################
#Apply the BUSgibbs Function
###############################################################################
set.seed(123)
BUSfits <- BUSgibbs(example_Data, n.subtypes = 3, n.iterations = 100, showIteration = FALSE)
summary(BUSfits)


\dontrun{
#estimated subtypes
est_subtypes <- Subtypes(BUSfits)

#estimated subtype effects
est_subtype_effects <- subtype_effects(BUSfits)

#estimated location batch effects
est_location_batch_effects <- location_batch_effects(BUSfits)

#BIC
BIC_val <- BIC(BUSfits)

#adjusted expression values
adjusted_data <- adjusted_values(BUSfits, Data)

#estimate intrinsic gene indicators by using the default postprob_DE_threshold = 0.5
est_L <- estimate_IG_indicators(BUSfits, postprob_DE_threshold = 0.5)

#obtain the intrinsic gene indicators
intrinsic_gene_indices <- IG_index(est_L)}

}
