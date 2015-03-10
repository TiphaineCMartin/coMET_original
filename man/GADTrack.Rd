\name{GADTrack}
\alias{GADTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one track of the genomic positions of variants from the Genetic Association Database (GAD) 
}
\description{
Create one track of the genomic positions of variants from the Genetic Association Database (GAD) (archive of human genetic association studies of complex diseases and disorders) using the Gviz bioconductor package
}
\usage{
GADTrack(gen, chr, start, end, showId=FALSE)
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
http://bioconductor.org/packages/release/bioc/html/Gviz.html

http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=gad
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
 \code{\link{genesNameENSEMBL}},  \code{\link{GeneReviewsTrack}}, \code{\link{genesENSEMBL}},  \code{\link{xenorefGenesUCSC}}, \code{\link{transcriptENSEMBL}},
}
\examples{
library("Gviz")
gen2 <- "hg19"
chrom2 <- "chr2"
start2 <- 38290160
end2 <- 38303219

if(interactive()) {
  gadtrack<-GADTrack(gen=gen2 ,chr=chrom2 ,start=start2 ,end=end2)
  plotTracks(gadtrack, from = start2, to =end2)
} else {
  data(gadtrack)
   plotTracks(gadtrack, from = start2, to =end2)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
