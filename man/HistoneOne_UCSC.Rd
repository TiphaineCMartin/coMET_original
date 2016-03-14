\name{HistoneOne_UCSC}
\alias{HistoneOne_UCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one track of one histone modification profile from the UCSC genome browser
}
\description{
Create one track of one histone modification profile from the UCSC genome browser (ENCODE/Broad) using the Gviz bioconductor package
}
\usage{
HistoneOne_UCSC(gen, chr, start, end, mySession, track.name = "Broad Histone", table.name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gen}{
the name of the genome. Data is not currently available for GRCh38 (hg38).
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
  \item{mySession}{
the object session from the function browserSession of rtracklayer
}
 \item{track.name}{
the name of the track, for example: "Broad Histone"
}
  \item{table.name}{
the name of the table from the track
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
An AnnotationTrack object of Gviz
}
\references{
http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=wgEncodeHistoneSuper

http://bioconductor.org/packages/release/bioc/html/Gviz.html
}
\author{
Tiphaine Martin
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{HistoneAll_UCSC}}
}
\examples{
library("Gviz")
library("rtracklayer")
gen <- "hg19"
chr <- "chr2"
start <- 38290160
end <- 38303219

if(interactive()) {
  BROWSER.SESSION="UCSC"
  mySession <- browserSession(BROWSER.SESSION)
  genome(mySession) <- gen
  histoneonetrack<-HistoneOne_UCSC(gen,chr,start,end,mySession)
  plotTracks(histoneonetrack, from = start, to =end)
} else {
  data(histoneonetrack)
  plotTracks(histoneonetrack, from = start, to =end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
