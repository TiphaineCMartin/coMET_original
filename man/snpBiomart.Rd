\name{snpBiomart}
\alias{snpBiomart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create track Short Variation from ENSEMBL
}
\description{
Create track Short Variation from ENSEMBL in using the Gviz bioconductor package
}
\usage{
snpBiomart(chr, start, end, dataset, showId=FALSE, title = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{chr}{
the chromosome of interest
}
  \item{start}{
the first position of region of interest (the smallest value)
}
  \item{end}{
the last position of region of interest (the biggest value)
} 
  \item{dataset}{
 The name of SNP database
}
  \item{showId}{
  Show the id of element or not
}
  \item{title}{
 The name of annotation track
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An AnnotationTrack object of Gviz
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tiphaine Martin
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library("Gviz")
gen <- "hg19"
chr <- "chr2"
start <- 38290160
end <- 38303219

if(interactive()){
  snptrack <- snpBiomart(chrom, start, end, 
                       dataset="hsapiens_snp_som",showId=FALSE)
  plotTracks(snptrack, from = start, to =end) 
} else {
  data(snpBiomarttrack)
  plotTracks(snptrack, from = start, to =end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
