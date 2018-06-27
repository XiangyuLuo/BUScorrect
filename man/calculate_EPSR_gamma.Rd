\name{calculate_EPSR_gamma}
\alias{calculate_EPSR_gamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the estimated potential scale reduction factors (EPSR factors) for location batch effects 
}
\description{
When Gibbs sampler attains stationary, the distances between multiple chains (with multiple initial values) should be small. The EPSR factors are calculated to help decide the iteration number of Gibbs sampler.
}
\usage{
calculate_EPSR_gamma(gamma_PosterSamp_chain1, gamma_PosterSamp_chain2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gamma_PosterSamp_chain1}{
	posterior samples of location batch effects from chain 1. 
}
  \item{gamma_PosterSamp_chain2}{
	posterior samples of location batch effects from chain 2.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{EPSR_gamma}{The EPSR factors for gamma, a G by B matrix. Note that EPSR_gamma[,1] is a NA vector. gamma[,1] are fixed at zero, so their EPSR factors are not taken into account. }
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
