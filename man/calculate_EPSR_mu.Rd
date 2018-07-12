\name{calculate_EPSR_mu}
\alias{calculate_EPSR_mu}
\title{
Calculate the estimated potential scale reduction factors (EPSR factors) for subtype effects
}
\description{
When Gibbs sampler attains stationary, the distances between multiple chains (with multiple initial values) should be small. The EPSR factors are calculated to help decide the iteration number of Gibbs sampler.
}
\usage{
calculate_EPSR_mu(mu_PosterSamp_chain1, mu_PosterSamp_chain2)
}

\arguments{
  \item{mu_PosterSamp_chain1}{
	posterior samples of subtype effects from chain 1. 
}
  \item{mu_PosterSamp_chain2}{
	posterior samples of subtype effects from chain 1. 
}
}

\value{
\item{EPSR_gamma}{The EPSR factors for mu, a G by K matrix. Note that EPSR_mu[,1] is a NA vector. mu[,1] are fixed at zero, so their EPSR factors are not taken into account. }

}

\author{
  Xiangyu Luo
}

\examples{
#10 genes, 2 subtypes, 100 posterior samples per parameter
chain1 <- 1+array(rnorm(10*2*100,sd=0.05), dim=c(10,2,100))
chain2 <- 1+array(rnorm(10*2*100,sd=0.05), dim=c(10,2,100))
calculate_EPSR_mu(chain1,chain2)
}

\keyword{ EPSR factors }
