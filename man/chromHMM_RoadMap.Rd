\name{chromHMM_RoadMap}
\alias{chromHMM_RoadMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a ChromHMM track from a file of RoadMap
}
\description{
Creates a ChromHMM track from a file of RoadMap using the Gviz bioconductor package. A complete list of features and their associated colours can be found in the user guide.
}
\usage{
chromHMM_RoadMap(gen="hg19",chr, start, end, bedFilePath, featureDisplay = 'all', colorcase='roadmap15')
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

  \item{featureDisplay}{
A vector of features to be displayed, such as 1_TssA. Spelling and capitalisation of features must be identical to those in the user guide (in the 'State & Acronym' column). There are three possibilities. First, the visualisation of only one feature (e.g. featureDisplay <- "1_TssA"), only the name of the specific feature is required. Second, visualisation of a set of features, for this a vector of features is required (e.g. featureDisplay <- c("1_TssA","2_TssAFlnk")). Finally, visualison all features in the genomic region, achived by using the word "all" (e.g. featureDisplay <- "all"), "all" is set by default. You can find the complete list of features and their associated colours in the user guide.
}
 \item{colorcase}{
the type of colors used to visualise different elements contained in ROADmap data with 15-,18-,25- states. choice between roadmap15, roadmap18, comet18, roadmap25 and comet25.


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
start <- 4500000
end <- 4600000
featureDisplay <- "7_Enh"

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
bedFilePath <- file.path(extdata, "RoadMap/E063_15_coreMarks_mnemonics.bed")

if(interactive()){
  chromHMM_RoadMapSingle <- chromHMM_RoadMap(gen="hg19",chr,start, end, bedFilePath, featureDisplay = featureDisplay, colorcase='roadmap15' )
  plotTracks(chromHMM_RoadMapSingle, from = start, to = end)
} else {
  data(chromHMM_RoadMapSingle)
  plotTracks(chromHMM_RoadMapSingle, from = start, to = end)
}

######

library("Gviz")
chr <- "chr1"
start <- 4500000
end <- 4600000
featureDisplay <- c("7_Enh","13_ReprPC")

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
bedFilePath <- file.path(extdata, "RoadMap/E063_15_coreMarks_mnemonics.bed")

if(interactive()){
  chromHMM_RoadMapMultiple <- chromHMM_RoadMap(gen="hg19",chr,start, end, bedFilePath, featureDisplay = featureDisplay, colorcase='roadmap15' )
  plotTracks(chromHMM_RoadMapMultiple, from = start, to = end)
} else {
  data(chromHMM_RoadMapMultiple)
  plotTracks(chromHMM_RoadMapMultiple, from = start, to = end)
}

#####

library("Gviz")
chr <- "chr1"
start <- 4500000
end <- 4600000
featureDisplay <- "all"

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
bedFilePath <- file.path(extdata, "RoadMap/E063_15_coreMarks_mnemonics.bed")

if(interactive()){
  chromHMM_RoadMapAll <- chromHMM_RoadMap(gen="hg19",chr,start, end, bedFilePath, featureDisplay = featureDisplay, colorcase='roadmap15' )
  plotTracks(chromHMM_RoadMapAll, from = start, to = end)
} else {
  data(chromHMM_RoadMapAll)
  plotTracks(chromHMM_RoadMapAll, from = start, to = end)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
