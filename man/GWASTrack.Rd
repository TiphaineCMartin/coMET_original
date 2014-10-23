\name{GWASTrack}
\alias{GWASTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one track of the genomic positions of variants from the GWAS catalog
}
\description{
Create one track of the genomic positions of variants from the NHGRI Catalog of Published Genome-Wide Association Studies using the Gviz bioconductor package
}
\usage{
GWASTrack(gen, chr, start, end, showId=FALSE)
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
 \item{showId}{
Show the ID of the genetic elements
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An UcscTrack object of Gviz
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
start <- 37949607
end <- 37965207

if(interactive()) {
  gwastrack <-GWASTrack(gen,chrom,start,end)
  plotTracks(gwastrack, from = start, to =end) 
} else {
  data(GWASTrack)
  plotTracks(gwastrack, from = start, to =end) 
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
