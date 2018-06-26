\name{summary.BUSfits}
\alias{summary.BUSfits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary the output by \code{BUSgibbs}
}
\description{Call the function to summarize the output object by \code{BUSgibbs}.}
\usage{
\method{summary}{BUSfits}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
The BUSfits object output by the function \code{BUSgibbs}.
}
  \item{\dots}{not used.}
}
\value{
	summarize the results from the output by BUSgibbs.
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

BUSfits <- BUSgibbs(example_Data, n.subtypes = 3, n.iterations = 100, showIteration = FALSE, seed=123)
summary(BUSfits)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

