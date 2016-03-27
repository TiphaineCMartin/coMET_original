\name{transcript_ENSEMBL}
\alias{transcript_ENSEMBL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a track of transcripts from ENSEMBL
}
\description{
Create a track to visualize different transcripts from ENSEMBL using the Gviz bioconductor package
}
\usage{
transcript_ENSEMBL(gen, chr, start, end,showId = FALSE)
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
 \code{\link{ISCA_UCSC}},  \code{\link{GWAScatalog_UCSC}}, \code{\link{knownGenes_UCSC}}, 
 \code{\link{genesName_ENSEMBL}},  \code{\link{GeneReviews_UCSC}}, \code{\link{GAD_UCSC}},  \code{\link{genes_ENSEMBL}}, \code{\link{xenorefGenes_UCSC}},
}
\examples{
library("Gviz")
gen <- "hg19"
chr <- "chr2"
start <- 32290160
end <- 33303219

if(interactive()){
  transENSMBLtrack<-transcript_ENSEMBL(gen,chr,start,end,showId=TRUE)
  plotTracks(transENSMBLtrack, from=start, to=end)
} else {
  data(transENSMBLtrack)
  plotTracks(transENSMBLtrack, from=start, to=end)
}
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
