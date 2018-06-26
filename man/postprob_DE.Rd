\name{postprob_DE}
\alias{postprob_DE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the posterior probability of being differentially expressed
}
\description{Calculate the posterior probability of being differentially expressed for genes in subtypes k (k>=2) compared to subtype 1.}
\usage{
postprob_DE(BUSfits)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BUSfits}{
	The BUSfits object output by the function BUSgibbs.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{postprob_DE_matr}{the matrix of posterior probabilities of being differentially expressed }
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

BUSfits <- BUSgibbs(example_Data, n.subtypes = 3, n.iterations = 100, showIteration = FALSE, seed=123)
postprob_DE(BUSfits)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

