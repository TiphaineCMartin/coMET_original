\name{eQTL}
\alias{eQTL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a track from a file for eQTL data
}
\description{
Creates a track from a BED file for eQTL data using the Gviz bioconductor package. A complete list of features and their associated colours can be found in the user guide.
}
\usage{
eQTL(gen,chr, start, end, bedFilePath, featureDisplay, showId=FALSE,type_stacking="squish",just_group="above" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{gen}{
the name of the genome.
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
A vector of eQTL features to be displayed, such as SNP. Spelling and capitalisation of features must be identical to those in the user guide. There are three possibilities. First, the visualisation of only one feature (e.g. featureDisplay <- "CpG"), only the name of the specific feature is required. Second, visualisation of a set of features, for this a vector of features is required (e.g. featureDisplay <- c("SNP","CpG")). Finally, visualison all features in the genomic region, achived by using the word "all" (e.g. featureDisplay <- "all"), "all" is set by default. You can find the complete list of features and their associated colours in the user guide.
}
  \item{showId}{
Allows to visualise the Id of eQTL group.
}
  \item{type_stacking}{
Object of class"character", the stacking type of overlapping items on the final plot.One in c(hide, dense, squish, pack,full). More information cf the option "stacking" in Gviz
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
chr <- "chr15"
start <- 74889136
end <-   75018200
featureDisplay <- "SNP"
gen="hg19"

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
bedFilePath <- file.path(extdata, "eQTL.bed")

if(interactive()){
  eQTLTrackSingle <- eQTL(gen,chr,start, end, bedFilePath, featureDisplay = featureDisplay )
  plotTracks(eQTLTrackSingle, from = start, to = end)
} else {
  data(eQTLTrackSingle)
  plotTracks(eQTLTrackSingle, from = start, to = end)
}

###########

library("Gviz")
chr <- "chr15"
start <- 74889136
end <-   75018200
featureDisplay <- c("SNP","mRNA_pheno")
gen="hg19"

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
bedFilePath <- file.path(extdata, "eQTL.bed")

if(interactive()){
  eQTLTrackMultiple <- eQTL(gen,chr,start, end, bedFilePath, featureDisplay = featureDisplay )
  plotTracks(eQTLTrackMultiple, from = start, to = end)
} else {
  data(eQTLTrackMultiple)
  plotTracks(eQTLTrackMultiple, from = start, to = end)
}

######

library("Gviz")
chr <- "chr15"
start <- 74889136
end <-   75018200
featureDisplay <- "all"
gen="hg19"

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
bedFilePath <- file.path(extdata, "eQTL.bed")

if(interactive()){
  eQTLTrackAll <- eQTL(gen,chr,start, end, bedFilePath, featureDisplay = featureDisplay )
  plotTracks(eQTLTrackAll, from = start, to = end)
} else {
  data(eQTLTrackAll)
  plotTracks(eQTLTrackAll, from = start, to = end)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
