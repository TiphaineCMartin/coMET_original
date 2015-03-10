\name{chromatinHMMAll}
\alias{chromatinHMMAll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating multiple chromHMM tracks from the UCSC genome browser 
}
\description{
Create multiple chromHMM Broad tracks by connecting to the UCSC genome browser using the GViz bioconductor package
}
\usage{
chromatinHMMAll(gen, chr, start, end, mySession, track.name = "Broad ChromHMM", 
                pattern = NULL, table.name = NULL)
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
the name of the track, for example : Broad ChromHMM
}
    \item{pattern}{
the pattern of the track to visualise
}
    \item{table.name}{
the name of the table from the track
}
 
}
%%\details{
%%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
    list of AnnotationTrack objects of GViz
}
\references{
http://bioconductor.org/packages/release/bioc/html/Gviz.html

http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=wgEncodeHistoneSuper
}
\author{
Tiphaine Martin
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{chromatinHMMOne}}
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
        track.name="Broad ChromHMM"
        tablestrack<-tableNames(ucscTableQuery(mySession, track=track.name))
        table.name<-tablestrack[1]
        PATTERN.REGULATION<-"GM12878"
        
        chromhmmPattern<-chromatinHMMAll(gen,chr,start,end,mySession,track.name,PATTERN.REGULATION)
        plotTracks(chromhmmPattern, from = start, to =end)

        chromhmmNoPattern<-chromatinHMMAll(gen,chr,start,end,mySession,track.name)
        plotTracks(chromhmmNoPattern, from = start, to =end)
    } else {

        data(chromhmmPattern)
        plotTracks(chromhmmPattern, from = start, to =end)

        data(chromhmmNoPattern)
        plotTracks(chromhmmNoPattern, from = start, to =end)
    }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

