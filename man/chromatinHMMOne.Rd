\name{chromatinHMMOne}
\alias{chromatinHMMOne}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating one chromHMM track from the UCSC genome browser
}
\description{
Create one track of only one type of chromHMM Broad element from the UCSC genome browser using the Gviz bioconductor package
}
\usage{
chromatinHMMOne(gen, chr, start, end, mySession, track.name = "Broad ChromHMM", 
                table.name = NULL)
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
the first position in region of interest (the smallest value)
}
  \item{end}{
the last position in region of interest (the biggest value)
}
  \item{mySession}{
the object session from the function browserSession of rtracklayer
}
  \item{track.name}{
the name of the track( Broad ChromHMM )
}
  \item{table.name}{
the name of the table from the track
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An AnnotationTrack object of Gviz
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
  end <- 38303219
if(interactive()) {

  BROWSER.SESSION="UCSC"
  mySession <- browserSession(BROWSER.SESSION)
  genome(mySession) <- gen
  track.name="Broad ChromHMM"
  tablestrack<-tableNames(ucscTableQuery(mySession, track=track.name))
  table.name<-tablestrack[1]
  chromhmmtrackone<-chromatinHMMOne(gen,chr,start,end,mySession,track.name,table.name)
  plotTracks(chromhmmtrackone, from = start, to =end) 
}else {

  data(chromhmmtrackone)
  plotTracks(chromhmmtrackone, from = start, to =end) 
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

