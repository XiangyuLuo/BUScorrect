\name{estimate_IG_indicators}
\alias{estimate_IG_indicators}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the intrinsic gene indicators
}
\description{Call the function to estimate the intrinsic gene indicators.}
\usage{
estimate_IG_indicators(BUSfits, postprob_DE_threshold = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BUSfits}{
The BUSfits object output by the function \code{BUSgibbs}.
}
  \item{postprob_DE_threshold}{
the threshold to call an intrinsic gene indicator to be one or not according to whether its posterior probability is higher than postprob_DE_threshold or not. The default is 0.5.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{est_L}{the estimated intrinsic gene indicators, a matrix where the rows represent genes and the columns correspond to subtypes k=2,...,K}
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

set.seed(123)
BUSfits <- BUSgibbs(example_Data, n.subtypes = 3, n.iterations = 100, showIteration = FALSE)
#select the posterior probability threshold to estimate the intrinsic gene indicators
thr0 <- postprob_DE_thr_fun(BUSfits, fdr_threshold=0.1)
est_L <- estimate_IG_indicators(BUSfits, postprob_DE_threshold=thr0)

#obtain the intrinsic gene indicators
intrinsic_gene_indices <- IG_index(est_L)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

