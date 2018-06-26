\name{calculate_EPSR_mu}
\alias{calculate_EPSR_mu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the estimated potential scale reduction factors (EPSR factors) for subtype effects
}
\description{
When Gibbs sampler attains stationary, the distances between multiple chains (with multiple initial values) should be small. The EPSR factors are calculated to help decide the iteration number of Gibbs sampler.
}
\usage{
calculate_EPSR_mu(mu_PosterSamp_chain1, mu_PosterSamp_chain2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu_PosterSamp_chain1}{
	posterior samples of subtype effects from chain 1. 
}
  \item{mu_PosterSamp_chain2}{
	posterior samples of subtype effects from chain 1. 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{EPSR_gamma}{The EPSR factors for mu, a G by K matrix. Note that EPSR_mu[,1] is a NA vector. mu[,1] are fixed at zero, so their EPSR factors are not taken into account. }

}
%\references{
%% ~put references to the literature/web site here ~
%}
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
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EPSR factors }
