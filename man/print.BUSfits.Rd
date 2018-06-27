\name{print.BUSfits}
\alias{print.BUSfits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print the output by BUSgibbs
}
\description{Call the function to print the output by BUSgibbs.}
\usage{
\method{print}{BUSfits}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The BUSfits object output by the function \code{BUSgibbs}.
}
  \item{\dots}{not used.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
	print the results from the output by BUSgibbs.
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
print(BUSfits)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
