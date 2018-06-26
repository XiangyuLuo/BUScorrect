\name{BUSgibbs}
\alias{BUSgibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Batch Effects Correction and Subtype Discovery using Gibbs Sampler 
}
\description{
The function "BUSgibbs" stands for fitting the Batch effects correction with Unknown Subtypes model (BUS) with the Gibbs Sampler. BUS is capable of (a) correcting batch effects explicitly, (b) grouping samples that share similar characteristics into subtypes, (c) identifying features that distinguish subtypes, and (d) enjoying a linear-order computation complexity. After correcting the batch effects with BUS, the corrected value can be used for other analysis as if all samples are measured in a single batch. We adopt the Bayesian framework and use the Gibbs sampler to conduct posterior inference for the BUS model. 
}
\usage{
BUSgibbs(Data, n.subtypes, n.iterations = 500, n.records = floor(n.iterations/2), 
		 hyperparameters = c(1, sqrt(5), sqrt(5), 2, 2, 1, 2, 0.005, 1, 3, 10), 
		 seed = 123, showIteration = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
	\code{Data} is an R list with the length equal to the batch number. The \code{b}th element of \code{Data} is the gene expression matrix in batch \code{b}, where the rows correspond to genes and the columns represent samples. 
}
  \item{n.subtypes}{
	\code{n.subtypes} is the subtype number, which needs to be specified by the user.
}
  \item{n.iterations}{
	\code{n.iterations} is the iteration number used in the Gibbs sampler. The default is 500.
}
  \item{n.records}{
	The posterior samples in the last \code{n.records} iterations are used to conduct posterior inference. The default is one half of \code{n.iterations}.
}
  \item{hyperparameters}{
	\code{hyperparameters} is a hyper-parameter vector with 11 elements used in the Gibbs sampler. The first element to the last element of \code{hyperparameters} are as follows.
	\code{eta_alpha}: the mean of the normal prior for alpha_g; 
	\code{tau_alpha}: the standard deviation of the normal prior for alpha_g;
	\code{tau_gamma}: the standard deviation of the normal prior for gamma;
	\code{alpha_par}: the parameter of the Dirichlet prior for subtype proportions;
	\code{a_inv_gamma}: the shape of the gamma prior for 1/sigma^2_{bg};
	\code{b_inv_gamma}: the rate of the gamma prior for 1/sigma^2_{bg};
	\code{a_tau0}: the shape of the gamma prior for 1/tau^2_{mu 0};
	\code{b_tau0}: the rate of the gamma prior for 1/tau^2_mu 0;
	\code{(a_p, b_p)}: parameters in the beta prior for p;
	\code{tau_{mu1}}: the standard deviation of the normal prior of mu_{gk} (k >= 2) when the gene expression level in subtype k is different from that in subtype one.
}
  \item{seed}{
	Different \code{seed} values allow the user to run the Gibbs sampler with different initial values.
}
  \item{showIteration}{
	If TRUE, the iteration number will be displayed when conducting Gibbs sampler. The default is TRUE. 
}
}
\details{
 Notice that \code{Data}, the input original gene expression values, are organized in the format of an R list with length equal to the batch number. Its \code{b}th element \code{Data[[b]]} is a \code{G} by \code{n_b} matrix, where \code{G} is the gene number and \code{n_b} is the sampler size of batch \code{b}. 
}
\value{
  \item{L_PosterSamp}{The posterior samples of the intrinsic gene indicators. The return is a G by K-1 by \code{n.records} array, where G is the gene number, K is the subtype number, and \code{n.records} is the number for recorded iterations.   }
  \item{Subtypes}{The estimated subtypes, an R list with length equal to the batch number B, in which Subtypes[[b]] is an integer vector showing the subtype indicators of samples in batch b.}
   \item{tau_mu_zero}{The estimated tau_{mu 0}, which is the prior normal distribution's standard deviation of the subtype effects when there is no differential expression.}
   \item{p}{The estimated proportion of intrinsic genes.}
   \item{pi}{The estimated subtype proportions across batches, a B by K matrix, whose [b,k] element is the estimated proportion of subtype k in the batch b.}
   \item{alpha}{The estimated baseline expression levels, a G-dimension vector, whose gth element is the estimated mean gene expression level of gene g in subtype one.}
   \item{gamma_PosterSamp}{The posterior samples of location batch effects, a G by B by \code{n.records} array.}
   \item{gamma}{The estimated location batch effects, a G by B matrix, where gamma_{gb} is the ``location'' batch effect on gene g in the batch b. Note that the first column is zero as the first batch is taken as the reference batch without batch effects.}
   \item{sigma_sq_PosterSamp}{The posterior samples of variances, a G by B by \code{n.records} array.}   
   \item{sigma_sq}{The estimated variance, a G by B matrix, whose [g,b] element is the variance of gene g's expression in the batch b.}
   \item{mu_PosterSamp}{The posterior samples of subtype effects, a G by K by \code{n.records} array.}
   \item{mu}{The estimated subtype effects, a G by K matrix, whose [g,k] element is the subtype k effect on gene g. Note that the first column is zero as the fist subtype is taken as the baseline subtype.}
   \item{BIC}{the BIC value when K = \code{n.subtypes}, which is used to determine the subtype number by varying the value of K.}
}
\references{
Xiangyu Luo, Yingying Wei. Batch Effects Correction with Unknown Subtypes. Journal of the American Statistical Association. Accepted.
}
\author{
  Xiangyu Luo
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
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
# pi[b,k] stands for the proportion of kth subtype in bth batch

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

BUSfits <- BUSgibbs(example_Data, n.subtypes = 3, n.iterations = 100, showIteration = FALSE, seed=123)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Batch Effects }
\keyword{ Subtype Discovery }
\keyword{ Gibbs Sampler}
\keyword{ Data Integration}
