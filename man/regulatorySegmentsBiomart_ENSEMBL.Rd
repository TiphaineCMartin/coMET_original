\name{regulatorySegmentsBiomart_ENSEMBL}
\alias{regulatorySegmentsBiomart_ENSEMBL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a binding motif track from ENSEMBL
}
\description{
Creates a track of regulatory segments from ENSEMBL using the Gviz bioconductor package. A complete list of features and their associated colours can be found in the user guide.
}
\usage{
regulatorySegmentsBiomart_ENSEMBL(gen, chr, start, end, featureDisplay = 'all', datasetEnsembl = "hsapiens_segmentation_feature")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{gen}{
The name of the genome. Currently only handles human data from either the previous version, GRCh37 (also known as hg19) or the current version, GRCh38 (also known as hg38).
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
  \item{featureDisplay}{
A vector of regulatory features to be displayed, such as Predicted heterochomatin. Spelling and capitalisation of features must be identical to those in the user guide. There are three possibilities. First, the visualisation of only one feature (e.g. featureDisplay <- "Predicted heterochomatin"), only the name of the specific feature is required. Second, visualisation of a set of features, for this a vector of features is required (e.g. featureDisplay <- c("Predicted low activity","Predicted heterochomatin")). Finally, visualison all features in the genomic region, achived by using the word "all" (e.g. featureDisplay <- "all"), "all" is set by default. You can find the complete list of features and their associated colours in the user guide.
}
  \item{datasetEnsembl}{
Allows the user to manually set which data set is used if required.Default=hsapiens_segmentation_feature
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
gen <- "hg38"
chr <- "chr1"
start <- 10000
end <- 50000
featureDisplay <- "CTCF enriched"

if(interactive()){
 regulatorySegmentsBiomartTrackSingle<-regulatorySegmentsBiomart_ENSEMBL(gen,chr,start,end,featureDisplay)
  plotTracks(regulatorySegmentsBiomartTrackSingle, from = start, to = end)
} else {
  data(regulatorySegmentsBiomartTrackSingle)
  plotTracks(regulatorySegmentsBiomartTrackSingle, from = start, to = end)
}

####

library("Gviz")
gen <- "hg38"
chr <- "chr1"
start <- 10000
end <- 50000
featureDisplay <- c("CTCF enriched","Predicted Promoter Flank")

if(interactive()){
 regulatorySegmentsBiomartTrackMultiple<-regulatorySegmentsBiomart_ENSEMBL(gen,chr,start,end,featureDisplay)
  plotTracks(regulatorySegmentsBiomartTrackMultiple, from = start, to = end)
} else {
  data(regulatorySegmentsBiomartTrackMultiple)
  plotTracks(regulatorySegmentsBiomartTrackMultiple, from = start, to = end)
}

####

library("Gviz")
gen <- "hg38"
chr <- "chr1"
start <- 10000
end <- 50000
featureDisplay <- "all"
if(interactive()){
 regulatorySegmentsBiomartTrackAll<-regulatorySegmentsBiomart_ENSEMBL(gen,chr,start,end,featureDisplay)
  plotTracks(regulatorySegmentsBiomartTrackAll, from = start, to = end)
} else {
  data(regulatorySegmentsBiomartTrackAll)
  plotTracks(regulatorySegmentsBiomartTrackAll, from = start, to = end)

}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
