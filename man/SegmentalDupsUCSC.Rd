\name{SegmentalDupsUCSC}
\alias{SegmentalDupsUCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one track of the genomic positions of regions from SegmentalDupsUCSC
}
\description{
Create one track of the genomic positions of regions from SegmentalDupsUCSC using the Gviz bioconductor package
}
\usage{
SegmentalDupsUCSC(gen, chr, start, end, showId=FALSE)
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

http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=rmsk
}
\author{
Tiphaine Martin

Tom Hardiman
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
library("Gviz")
gen <- "hg19"
chr <- "chr2"
start <- 100000
end <- 200000

DupTrack <-SegmentalDupsUCSC(gen,chr,start,end, showID=TRUE)
plotTracks(DupTrack, from = start, to = end)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
