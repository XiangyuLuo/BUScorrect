\name{visualize_data}
\alias{visualize_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize the gene expression data from multiple batches
}
\description{Use "heatmap.2" in R package "gplots" to visualize the gene expression data across multiple batches.}
\usage{
visualize_data(Data, title_name="Heatmap", gene_ind_set, color_key_range=seq(-0.5,8.5,1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
	The gene expression data, an R list with length equal to the batch number. Each of its element is a gene expression matrix, where rows are genes and columns represent samples.
	}
   \item{title_name}{
	The title name of the heatmap.
   }
   \item{gene_ind_set}{The indices of the set of genes the user wants to display in the heatmap.}
   \item{color_key_range}{The color range in the color key.}
}
\details{
  The values displayed in the heatmap are the raw values in the argument \code{Data} without scaling.
}
\value{
	visualize the gene expression data matrix, where one row is a gene and one column represents a sample.
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

visualize_data(example_Data, title_name="Heatmap", gene_ind_set = 1:20, color_key_range=seq(0,10,2))

}

\keyword{ visualization }% __ONLY ONE__ keyword per line
