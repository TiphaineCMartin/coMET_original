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
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
An UcscTrack object of Gviz
}
\references{
http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=gwasCatalog

http://bioconductor.org/packages/release/bioc/html/Gviz.html
}
\author{
Tiphaine Martin
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{ISCATrack}},  \code{\link{knownGenesUCSC}}, \code{\link{genesNameENSEMBL}}, 
 \code{\link{GeneReviewsTrack}},  \code{\link{GADTrack}}, \code{\link{genesENSEMBL}},  \code{\link{xenorefGenesUCSC}}, \code{\link{transcriptENSEMBL}},
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
