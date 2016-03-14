\name{gcContent_UCSC}
\alias{gcContent_UCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one track of GC content from UCSC
}
\description{
Create a track of GC content from UCSC using the Gviz bioconductor package
}
\usage{
gcContent_UCSC(gen, chr, start, end)
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
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
A UcscTrack object of Gviz
}
\references{
http://bioconductor.org/packages/release/bioc/html/Gviz.html

http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=gc5Base
}
\author{
Tiphaine Martin
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
library("Gviz")
gen <- "hg38"
chr <- "chr7"
start <- 38290160
end <- 38303219

if(interactive()){
  gctrack<-gcContent_UCSC(gen,chr,start,end)
  plotTracks(gctrack,from= start, to=end)
} else {
  data(gctrack)
  plotTracks(gctrack,from= start, to=end)
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

