\name{calculate_EPSR_sigma_sq}
\alias{calculate_EPSR_sigma_sq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the estimated potential scale reduction factors (EPSR factors) for the variances of expression values 
}
\description{
When Gibbs sampler attains stationary, the distances between multiple chains (with multiple initial values) should be small. The EPSR factors are calculated to help decide the iteration number of Gibbs sampler.
}
\usage{
calculate_EPSR_sigma_sq(sigma_sq_PosterSamp_chain1, sigma_sq_PosterSamp_chain2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigma_sq_PosterSamp_chain1}{
	posterior samples of variances from chain 1. 
}
  \item{sigma_sq_PosterSamp_chain2}{
	posterior samples of variances from chain 2.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{EPSR_sigma_sq}{The EPSR factors for sigma_sq, a G by B matrix.}
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
