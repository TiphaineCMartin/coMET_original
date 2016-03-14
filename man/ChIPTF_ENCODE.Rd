\name{ChIPTF_ENCODE}
\alias{ChIPTF_ENCODE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a TF motif track from ENCODE
}
\description{
Creates a track of TF motifs from ENCODE using the Gviz bioconductor package. A complete list of features and their associated colours can be found in the user guide.
}
\usage{
ChIPTF_ENCODE(gen="hg19", chr, start, end, bedFilePath, featureDisplay='all', motifColorFile, type_stacking='dense',showId=FALSE,just_group="above")
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
The path of the BED file from Kheradpour and Kellis, 2014.
}
  \item{featureDisplay}{
A vector of regulatory features to be displayed, such as Predicted heterochomatin. Spelling and capitalisation of features must be identical to those in the user guide. There are three possibilities. First, the visualisation of only one feature (e.g. featureDisplay <- "Predicted heterochomatin"), only the name of the specific feature is required. Second, visualisation of a set of features, for this a vector of features is required (e.g. featureDisplay <- c("Predicted low activity","Predicted heterochomatin")). Finally, visualison all features in the genomic region, achived by using the word "all" (e.g. featureDisplay <- "all"), "all" is set by default. You can find the complete list of features and their associated colours in the user guide.
}
\item{motifColorFile}{
The path of the BED file with 2 columns ( the first for motif name and the second for the color in hex format without \# in the beginning) with a header.
}
  \item{type_stacking}{
Object of class"character", the stacking type of overlapping items on the final plot.One in c(hide, dense, squish, pack,full). More information cf the option "stacking" in Gviz
}
\item{showId}{
logical. say if we write the name of group
}
\item{just_group}{
position. say where we write the name of group (choice in c("above","righ","left"))
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

Got to BindingMotifsBiomart binding motif biomart
}
\author{
Tiphaine Martin
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
gen <- "hg19"
chr<-"chr1"
start <- 1000
end <- 329000


if(interactive()){
    extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
    bedFilePath <- file.path(extdata, "ENCODE/motifs1000_matches_ENCODE.txt")
    motif_color <- file.path(extdata, "ENCODE/TFmotifs_colors.csv")
    chipTFtrack <- ChIPTF_ENCODE(gen,chr,start, end, bedFilePath, featureDisplay=c("AHR::ARNT::HIF1A_1","AIRE_1","AIRE_2","AHR::ARNT_1"), motif_color,type_stacking="squish",showId=TRUE)
    plotTracks(chipTFtrack, from = start, to = end)
} else {
    data(chipTFtrack)
    plotTracks(chipTFtrack, from = start, to = end)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
