\name{regulationBiomart}
\alias{regulationBiomart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a regulation track from ENSEMBL
}
\description{
Create a 'Regulation' track from ENSEMBL using the Gviz bioconductor package
}
\usage{
regulationBiomart(gen, chr, start, end)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{gen}{
the name of the genome
}
  \item{chr}{
the chromosome of interest
}
  \item{start}{
the first position in the region of interest (the smallest value)
}
  \item{end}{
the last position in the region of interest (the largest value)
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
chr <- "chr7"
start <- 38290160
end <- 38303219

if(interactive()){
  regulationENSEMBLtrack<-regulationBiomart(gen,chr,start,end)
  plotTracks(regulationENSEMBLtrack, from = start, to =end) 
} else {
  data(regulationENSEMBLtrack)
  plotTracks(regulationENSEMBLtrack, from = start, to =end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
