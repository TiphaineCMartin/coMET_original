\name{metQTL}
\alias{metQTL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a track from a file for metQTL data
}
\description{
Creates a track from a BED file for metQTL data using the Gviz bioconductor package. A complete list of features and their associated colours can be found in the user guide.
}
\usage{
metQTL(gen, chr, start, end, bedFilePath, featureDisplay, showId=FALSE,type_stacking="squish",just_group="above")
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
A vector of metQTL features to be displayed, such as SNP. Spelling and capitalisation of features must be identical to those in the user guide. There are three possibilities. First, the visualisation of only one feature (e.g. featureDisplay <- "CpG"), only the name of the specific feature is required. Second, visualisation of a set of features, for this a vector of features is required (e.g. featureDisplay <- c("SNP","CpG")). Finally, visualison all features in the genomic region, achived by using the word "all" (e.g. featureDisplay <- "all"), "all" is set by default. You can find the complete list of features and their associated colours in the user guide.
}
  \item{showId}{
Allows the visualization of the Id of metQTL group.
}
  \item{type_stacking}{
Sets the type of stacking used by Gviz for plots. By default this is set to 'squish'. For more information see Gviz user guide.
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

Got to ENSEMBLregulation binding motif biomart
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

gen <- 'hg19'
chr <- "chr15"
start <- 74889136
end <-   75018200
featureDisplay <- "trans_local_metQTL"
type_stacking <- "squish"

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
mqtlbedFilePath <- file.path(extdata, "metQTL.bed")

if(interactive()){
  metQTLTrackSingle <- metQTL(gen,chr,start, end,mqtlbedFilePath, featureDisplay = featureDisplay )
  plotTracks(metQTLTrackSingle, from = start, to = end)
} else {
  data(metQTLTrackSingle)
  plotTracks(metQTLTrackSingle, from = start, to = end)
}

###

library("Gviz")

gen <- 'hg19'
chr <- "chr15"
start <- 74889136
end <-   75018200

featureDisplay <- c("trans_local_metQTL","CpG")

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
bedFilePath <- file.path(extdata, "metQTL.bed")

if(interactive()){
  metQTLTrackMultiple <- metQTL(gen,chr,start, end, bedFilePath, featureDisplay = featureDisplay )
  plotTracks(metQTLTrackMultiple, from = start, to = end)
} else {
  data(metQTLTrackMultiple)
  plotTracks(metQTLTrackMultiple, from = start, to = end)
}

#####

library("Gviz")

gen <- 'hg19'
chr <- "chr15"
start <- 74889136
end <-   75018200

featureDisplay <- "all"

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
bedFilePath <- file.path(extdata, "metQTL.bed")

if(interactive()){
  metQTLTrackAll <- metQTL(gen,chr,start, end, bedFilePath, featureDisplay = featureDisplay )
  plotTracks(metQTLTrackAll, from = start, to = end)
} else {
  data(metQTLTrackAll)
  plotTracks(metQTLTrackAll, from = start, to = end)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
