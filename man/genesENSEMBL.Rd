\name{genesENSEMBL}
\alias{genesENSEMBL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one track of the genes in the genomic regions of interest from EMSEMBL
}
\description{
Create one track of the genes in the genomic regions of interest from EMSEMBL using the Gviz bioconductor package
}
\usage{
genesENSEMBL(gen, chr, start, end, showId=FALSE)
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
A BiomartGeneRegionTrack object of Gviz
}
\references{
http://bioconductor.org/packages/release/bioc/html/Gviz.html

http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=ensGene
}
\author{
Tiphaine Martin
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{ISCATrack}},  \code{\link{GWASTrack}}, \code{\link{knownGenesUCSC}}, 
 \code{\link{genesNameENSEMBL}},  \code{\link{GeneReviewsTrack}}, \code{\link{GADTrack}},  \code{\link{xenorefGenesUCSC}}, \code{\link{transcriptENSEMBL}},
}
\examples{
library("Gviz")
gen <- "hg19"
chrom <- "chr2"
start <- 38290160
end <- 38303219
if(interactive()) {
  genetrack <-genesENSEMBL(gen,chrom,start,end,showId=TRUE)
  plotTracks(genetrack, from = start, to =end) 
} else {
  data(geneENSEMBLtrack)
  plotTracks(genetrack, from = start, to =end)
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

