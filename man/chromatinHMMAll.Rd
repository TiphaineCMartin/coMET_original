\name{chromatinHMMAll}
\alias{chromatinHMMAll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create multiple chromaHMM Broad tracks from UCSC's genome browser 
}
\description{
Create multiple chromaHMM Broad tracks from a connection to UCSC genome browser in using the GViz bioconductor package
}
\usage{
chromatinHMMAll(gen, chr, start, end, mySession, track.name = NULL, pattern = NULL,  table.name = NULL)
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
   \item{track.name}{
the name of track, for example : Broad ChromHMM
}
  \item{pattern}{
the pattern of track to visualise
}
  \item{table.name}{
the name of table from the track
}
 
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 list of AnnotationTrack objects of GViz
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
require("Gviz")
gen <- "hg19"
BROWSER.SESSION="UCSC"
mySession <- browserSession(BROWSER.SESSION)
genome(mySession) <- gen
track.name="Broad ChromHMM"
tablestrack<-tableNames(ucscTableQuery(mySession, track=track.name))
table.name<-tablestrack[1]
chr <- "chr7"
start <- 38290160
end <- 38313219
PATTERN.REGULATION<-"GM12878"
tmp<-chromatinHMMAll(gen,chr,start,end,mySession,track.name,PATTERN.REGULATION)
plotTracks(tmp)

tmp<-chromatinHMMAll(gen,chr,start,end,mySession,track.name)
plotTracks(tmp, from = start, to =end)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

