\name{location_batch_effects}
\alias{location_batch_effects}

\title{
Obtain the location batch effects from the output by BUSgibbs
}
\description{Call the function to obtain the location batch effects from the output by \code{BUSgibbs}.}
\usage{
location_batch_effects(BUSfits)
}

\arguments{
  \item{BUSfits}{
The BUSfits object output by the function \code{BUSgibbs}.
}
}

\value{
\item{est_location_batch_effects}{The estimated location batch effects, a matrix whose rows represent genes and columns correspond to batches. }
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
est_location_batch_effects <- location_batch_effects(BUSfits)

}

