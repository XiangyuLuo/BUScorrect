\name{postprob_DE_thr_fun}
\alias{postprob_DE_thr_fun}

\title{
Select the the posterior probability threshold to control the false discovery rate
}
\description{To control the false discovery rate at the targeted level, call postprob_DE_thr_fun to obtain the threshold for the posterior probability of being differentially expressed.}
\usage{
postprob_DE_thr_fun(BUSfits, fdr_threshold = 0.1)
}

\arguments{
  \item{BUSfits}{
The BUSfits object output by the function \code{BUSgibbs}.
}
  \item{fdr_threshold}{
the false discovery rate level we want to control.
}
}

\value{
\item{thre0}{the posterior probability threshold that controls the false discovery rate.}
}
\references{
Xiangyu Luo, Yingying Wei. Batch Effects Correction with Unknown Subtypes. Journal of the American Statistical Association. Accepted.
}
\author{
  Xiangyu Luo
}

\examples{
rm(list = ls(all = TRUE))  
set.seed(123)
#a toy example, there are 6 samples and 20 genes in each batch
example_Data <- list()

#batch 1
example_Data[[1]] <- rbind(matrix(c(1,1,5,5,10,10,
						3,3,7,7,12,12), ncol=6, byrow=TRUE), matrix(c(1,2),nrow=18, ncol=6))

#batch 2
batch2_effect <- c(2,2,2,1,1)
example_Data[[2]] <- rbind(matrix(c(1,1,5,5,10,10,
						3,3,7,7,12,12), ncol=6, byrow=TRUE), matrix(c(1,2),nrow=18, ncol=6)) + batch2_effect

#batch 3
batch3_effect <- c(3,2,1,1,2)
example_Data[[3]] <- rbind(matrix(c(1,1,5,5,10,10,
						3,3,7,7,12,12), ncol=6, byrow=TRUE), matrix(c(1,2),nrow=18, ncol=6)) + batch3_effect
set.seed(123)
BUSfits <- BUSgibbs(example_Data, n.subtypes = 3, n.iterations = 100, showIteration = FALSE)
#select kappa to estimate intrinsic gene indicators
thr0 <- postprob_DE_thr_fun(BUSfits, fdr_threshold=0.1)
est_L <- estimate_IG_indicators(BUSfits, postprob_DE_threshold=thr0)

#obtain the intrinsic gene indicators
intrinsic_gene_indices <- IG_index(est_L)
}


