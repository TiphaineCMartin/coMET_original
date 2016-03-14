\name{HiCdata2matrix}
\alias{HiCdata2matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a HiC matrix from a file (Rao et al., 2014)
}
\description{
Creates a HiC matrix from Rao et al.,2014.
}
\usage{
HiCdata2matrix( chr, start, end, bedFilePath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{
The chromosome of interest
}
  \item{start}{
The starting position in the region of interest (the smallest value)
}
  \item{end}{
The end position in the region of interest (the largest value)
}
 \item{bedFilePath}{
The path of the BED file from Kheradpour and Kellis, 2014.
}
}
%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
\value{
An AnnotationTrack object of Gviz
}
\references{
http://bioconductor.org/packages/release/bioc/html/Gviz.html

Got to BindingMotifsBiomart binding motif biomart
}
\author{
Tiphaine Martin
}
%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
library("corrplot")
gen <- "hg19"
chr<-"chr1"
start <- 5000000
end <- 9000000

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
bedFilePath <- file.path(extdata, "HiC/chr1_1mb.RAWobserved")

if(interactive()){
  matrix_HiC_Rao <- HiCdata2matrix(chr,start, end, bedFilePath)
cor_matrix_HiC <- cor(matrix_HiC_Rao)
diag(cor_matrix_HiC)<-1
corrplot(cor_matrix_HiC, method = "circle")
} else {
  data(matrix_HiC_Rao)
  cor_matrix_HiC <- cor(matrix_HiC_Rao)
diag(cor_matrix_HiC)<-1
corrplot(cor_matrix_HiC, method = "circle")
}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
