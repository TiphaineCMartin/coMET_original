\name{comet.list}
\alias{comet.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List the correlations between omic features
}
\description{
coMET is an R-based package to visualize EWAS (epigenome-wide association scans) results in a genomic region of interest. The main feature of coMET is to plot the the significance level of EWAS results in the selected region, along with correlation in DNA methylation values between CpG sites in the region. The coMET package generates plots of phenotype-association, co-methylation patterns, and a series of annotation tracks.
In addition, the function comet.list gives the list of correlations between omic features
}
\usage{
comet.list(cormatrix.file = NULL, cormatrix.method = "spearman", cormatrix.format = "raw", 
    cormatrix.conf.level=0.05, cormatrix.sig.level= 1, cormatrix.adjust="none", 
    cormatrix.type = "listdataframe", cormatrix.output="cormatrix_list",
    config.file = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{cormatrix.file}{
    Name of the raw data file or the pre-computed correlation matrix file. It is mandatory and has to be a file in tabular format with an header.
}
    \item{cormatrix.method}{
    Options for calculating the correlation matrix: spearman, pearson and kendall. Default value= spearman
}
    \item{cormatrix.format}{
    Format of the input cormatrix.file. TThere are two options: raw file (raw if CpG sites are by column and samples by row or raw_rev if CpG site are by row and samples by column) and pre-computed correlation matrix (cormatrix)
}
    \item{cormatrix.conf.level}{
    Alpha level for the confidence interval. Default value= 0.05. CI will be the alpha/2 lower and upper values.
}
    \item{cormatrix.sig.level}{
    Significant level to visualise the correlation. If the correlation has a pvalue below the significant level, the correlation will be colored in "goshwhite", else the color is related to the correlation level and the color scheme choosen.Default value =1.
}
    \item{cormatrix.adjust}{
    indicates which adjustment for multiple tests should be used. "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".Default value="none"
}
    \item{cormatrix.type}{
    Format of cormatrix.file. There are 2 different options: listfile or listdataframe.
}
    \item{cormatrix.output}{
    The path and the name of the output file without the extension
}
    \item{config.file}{
    Configuration file contains the values of these options instead of defining these by command line. It is a file where each line is one option. The name of option and its value are separated by "=".
}
    \item{verbose}{
       logical option TRUE or FALSE. TRUE (default). If TRUE, shows comments.
}
}

%\details{
%}
\value{
Create a list of correlation between omic features
}
\references{
http://epigen.kcl.ac.uk/comet/
}
\author{
Tiphaine Martin
}
%\note{
% the correlation test is done with the function cor.test of the package psych
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{comet.web}},\code{\link{comet}}
}
\examples{
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
myoutput <- file.path(extdata, "cyp1b1_res37_cormatrix_list_BH05.txt")


comet.list(cormatrix.file=mycorrelation,cormatrix.method = "spearman",
           cormatrix.format= "raw", cormatrix.conf.level=0.05,
           cormatrix.sig.level= 0.05, cormatrix.adjust="BH",
           cormatrix.type = "listfile", cormatrix.output=myoutput, 
           verbose=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% __ONLY ONE__ keyword per line
