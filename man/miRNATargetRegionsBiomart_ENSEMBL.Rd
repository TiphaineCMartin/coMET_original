\name{miRNATargetRegionsBiomart_ENSEMBL}
\alias{miRNATargetRegionsBiomart_ENSEMBL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a track of miRNA target regions from ENSEMBL
}
\description{
Creates a track of miRNA target regions from ENSEMBL using the Gviz bioconductor package.
}
\usage{
miRNATargetRegionsBiomart_ENSEMBL(gen, chr, start, end, showId=FALSE, datasetEnsembl = "hsapiens_mirna_target_feature")
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
 \item{showId}{
Show the ID of the genetic elements
}
  \item{datasetEnsembl}{
Allows the user to manually set which data set is used if required.Default=hsapiens_mirna_target_feature
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
start <- 1000000
end <- 20000000


if(interactive()){
 miRNATargetRegionsBiomartTrack<-miRNATargetRegionsBiomart_ENSEMBL(gen,chr,start,end,
                                      datasetEnsembl = "hsapiens_mirna_target_feature")
  plotTracks(miRNATargetRegionsBiomartTrack, from = start, to = end)
} else {
  data(miRNATargetRegionsBiomartTrack)
  plotTracks(miRNATargetRegionsBiomartTrack, from = start, to = end)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
