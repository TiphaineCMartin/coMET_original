\name{HistoneAll}
\alias{HistoneAll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create multiple tracks of Histone from UCSC's genome browser
}
\description{
Create multiple tracks of Histone from UCSC's genome browser (ENCODE/Broad) in using the Gviz bioconductor package
}
\usage{
HistoneAll(gen, chr, start, end, mySession, pattern = NULL, track.name = NULL, table.name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{gen}{
the name of genome
}
  \item{chr}{
the chromosome of interest
}
  \item{start}{
the first position of region of interest (the smallest value)
}
  \item{end}{
the last position of region of interest (the biggest value)
}
  \item{mySession}{
the object session from the function browserSession of rtracklayer
}
  \item{pattern}{
    The pattern of cell type
}
 \item{track.name}{
the name of track, for example : Histone
}
  \item{table.name}{
the name of table from the track
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of AnnotationTrack object of Gviz
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
end <- 38313219

if(interactive()){
  BROWSER.SESSION="UCSC"
  mySession <- browserSession(BROWSER.SESSION)
  genome(mySession) <- gen
  pattern1 <- "GM12878"

  histonalltrack<-HistoneAll(gen,chr,start,end,mySession, pattern=pattern1,track.name="Broad Histone")
  plotTracks(histonalltrack, from = start, to =end)
} else {
  data(histonalltrack)
  plotTracks(histonalltrack, from = start, to =end)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
