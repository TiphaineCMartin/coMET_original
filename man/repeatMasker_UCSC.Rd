\name{repeatMasker_UCSC}
\alias{repeatMasker_UCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one track of the genomic positions of regions from repeatMasker_UCSC
}
\description{
Create one track of the genomic positions of regions from repeatMasker_UCSC using the Gviz bioconductor package
}
\usage{
repeatMasker_UCSC(gen, chr, start, end, showId=FALSE,type_stacking="full")
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
 \item{type_stacking}{
the type of stacking data for this track. More information go to Gviz (the option "stacking")
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
gen <- "hg38"
chr <- "chr2"
start <- 38290160
end <- 38303219

if(interactive()){
  rmtrack <-repeatMasker_UCSC(gen,chr,start,end,showId=TRUE)
  plotTracks(rmtrack, from = start, to = end)
} else {
  data(repeatMaskerTrack)
  plotTracks(rmtrack, from = start, to = end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
