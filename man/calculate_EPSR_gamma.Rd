\name{calculate_EPSR_gamma}
\alias{calculate_EPSR_gamma}

\title{
Calculate the estimated potential scale reduction factors (EPSR factors) for location batch effects 
}
\description{
When Gibbs sampler attains stationary, the distances between multiple chains (with multiple initial values) should be small. The EPSR factors are calculated to help decide the iteration number of Gibbs sampler.
}
\usage{
calculate_EPSR_gamma(gamma_PosterSamp_chain1, gamma_PosterSamp_chain2)
}

\arguments{
  \item{gamma_PosterSamp_chain1}{
	posterior samples of location batch effects from chain 1. 
}
  \item{gamma_PosterSamp_chain2}{
	posterior samples of location batch effects from chain 2.
}
}

\value{
\item{EPSR_gamma}{The EPSR factors for gamma, a G by B matrix. Note that EPSR_gamma[,1] is a NA vector. gamma[,1] are fixed at zero, so their EPSR factors are not taken into account. }
}

\author{
  Xiangyu Luo
}


\examples{
#2 batches, 10 genes, 100 posterior samples per parameter
chain1 <- 1+array(rnorm(10*2*100,sd=0.05), dim=c(2,10,100))
chain2 <- 1+array(rnorm(10*2*100,sd=0.05), dim=c(2,10,100))
calculate_EPSR_gamma(chain1,chain2)
}

\keyword{ EPSR factors }
