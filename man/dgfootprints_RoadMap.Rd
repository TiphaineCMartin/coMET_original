\name{dgfootprints_RoadMap}
\alias{dgfootprints_RoadMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a track of DNA motif positional bias in digital genomic Footprinting Sites (DGFP) from a file of RoadMap
}
\description{
Creates a DGFP track from a file of RoadMap using the Gviz bioconductor package. A complete list of features and their associated colours can be found in the user guide.
}
\usage{
dgfootprints_RoadMap(gen="hg19", chr, start, end, bedFilePath, tissueGroupDisplay='Blood & T-cell',showId=FALSE, type_stacking="dense")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{gen}{
the name of the genome. Default value=hg19
}
  \item{chr}{
The chromosome of interest
}
  \item{start}{
The starting position in the region of interest (the smallest value)
}
  \item{end}{
The end position in the region of interest (the largest value)
}

  \item{bedFilePath}{
The file path to the .BED file containing the data to be visualised
}

  \item{tissueGroupDisplay}{
the group of tissue visualised among list("Neurosph","Epithelial","IMR90","Thymus","Heart","Brain","Digestive","Muscle","Other","iPSC","HSC & B-cell","Blood & T-cell"="ES-deriv")
}
\item{showId}{
logical. say if we write the name of group
}
  \item{type_stacking}{
Object of class"character", the stacking type of overlapping items on the final plot.One in c(hide, dense, squish, pack,full). More information cf the option "stacking" in Gviz
}
}
%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
\value{
An AnnotationTrack object of Gviz
}
\references{
http://bioconductor.org/packages/release/bioc/html/Gviz.html

Got to RoadMap Epigenome
}
\author{
Tiphaine Martin

Tom Hardiman
}
%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
library("Gviz")
chr <- "chr1"
start <- 236728
end <- 238778
gen="hg19"

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
bedFilePath <- file.path(extdata, "RoadMap/CD3-DS17198.hg19.bed")

if(interactive()){
  dgfootprints_RoadMapSingle <- dgfootprints_RoadMap(gen,chr,start, end, bedFilePath, tissueGroupDisplay='Blood & T-cell' )
  plotTracks(dgfootprints_RoadMapSingle, from = start, to = end)
} else {
  data(dgfootprints_RoadMapSingle)
  plotTracks(dgfootprints_RoadMapSingle, from = start, to = end)
}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
