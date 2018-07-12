\name{calculate_EPSR_sigma_sq}
\alias{calculate_EPSR_sigma_sq}
\title{
Calculate the estimated potential scale reduction factors (EPSR factors) for the variances of expression values 
}
\description{
When Gibbs sampler attains stationary, the distances between multiple chains (with multiple initial values) should be small. The EPSR factors are calculated to help decide the iteration number of Gibbs sampler.
}
\usage{
calculate_EPSR_sigma_sq(sigma_sq_PosterSamp_chain1, sigma_sq_PosterSamp_chain2)
}

\arguments{
  \item{sigma_sq_PosterSamp_chain1}{
	posterior samples of variances from chain 1. 
}
  \item{sigma_sq_PosterSamp_chain2}{
	posterior samples of variances from chain 2.
}
}

\value{
\item{EPSR_sigma_sq}{The EPSR factors for sigma_sq, a G by B matrix.}
}

\author{
  Xiangyu Luo
}

\examples{
#2 batches, 10 genes, 100 posterior samples per parameter
chain1 <- 1+array(rnorm(10*2*100,sd=0.05), dim=c(2,10,100))
chain2 <- 1+array(rnorm(10*2*100,sd=0.05), dim=c(2,10,100))
calculate_EPSR_sigma_sq(chain1,chain2)
}

\keyword{ EPSR factors }
