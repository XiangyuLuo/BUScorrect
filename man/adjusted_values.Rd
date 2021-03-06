\name{adjusted_values}
\alias{adjusted_values}

\title{
Obtain adjusted expression values from the output by BUSgibbs
}
\description{Call the function to obtain the corrected expression values with batch effects removed from the original input data.}
\usage{
adjusted_values(BUSfits, original_data)
}

\arguments{
  \item{BUSfits}{
The BUSfits object output by the function \code{BUSgibbs}.
}
  \item{original_data}{
The original gene expression data list with length equal to the batch number. Each of its element is a gene expression matrix for a specific batch, in which each row corresponds to a gene and each column corresponds to a sample.
}
}

\value{
\item{adjusted_data}{An R list with length equal to the batch number. Each of its element is a corrected gene expression matrix for a specific batch, in which each row corresponds to a gene and each column corresponds to a sample.}
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
adjusted_data <- adjusted_values(BUSfits, example_Data)

}
