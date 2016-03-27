\name{chromatinHMMOne_UCSC}
\alias{chromatinHMMOne_UCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating one chromHMM track from the UCSC genome browser
}
\description{
Create one track of only one type of chromHMM Broad element from the UCSC genome browser using the Gviz bioconductor package
}
\usage{
chromatinHMMOne_UCSC(gen, chr, start, end, mySession, color="coMET", table.name = NULL)
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
the first position in region of interest (the smallest value)
}
    \item{end}{
the last position in region of interest (the biggest value)
}
    \item{mySession}{
the object session from the function browserSession of rtracklayer
}
    \item{color}{
the color scheme used for plots. By defult this is set to 'coMET' to allow easy indentifcation of differnent elements. The color scheme set by UCSC can also be used. Consult userguide for table of colors.
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
http://bioconductor.org/packages/release/bioc/html/Gviz.html

http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=wgEncodeHistoneSuper
}
\author{
Tiphaine Martin

Tom Hardiman
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{chromatinHMMAll_UCSC}}
}
\examples{
    library("Gviz")
    library("rtracklayer")
    gen <- "hg19"
    chr <- "chr2"
    start <- 38290160
    end <- 38303219
    color <- "coMET"

    if(interactive()) {
        BROWSER.SESSION="UCSC"
        mySession <- browserSession(BROWSER.SESSION)
        genome(mySession) <- gen
        track.name="Broad ChromHMM"
        tablestrack<-tableNames(ucscTableQuery(mySession, track=track.name))
        table.name<-tablestrack[1]
        chromhmmtrackone<-chromatinHMMOne_UCSC(gen,chr,start,end,mySession,color="coMET",table.name)
        plotTracks(chromhmmtrackone, from = start, to =end)
    }else {
        data(chromhmmtrackone)
        plotTracks(chromhmmtrackone, from = start, to =end)
    }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

