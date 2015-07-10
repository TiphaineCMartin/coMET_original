\name{BindingMotifsBiomart}
\alias{BindingMotifsBiomart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a binding motif track from ENSEMBL
}
\description{
Creates a binding motif track from ENSEMBL using the Gviz bioconductor package
}
\usage{
BindingMotifsBiomart(gen, chr, start, end, featureDisplay = "all", datasetEnsembl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{gen}{
The name of the genome. Currently only manages human data from the previous version GRCh37, also called hg19, and the latest version GRCh38, also called hg38.
}
  \item{chr}{
The chromosome of interest
}
  \item{start}{
The first position in the region of interest (the smallest value)
}
  \item{end}{
The last position in the region of interest (the largest value)
}
  \item{featureDisplay}{
A vector of regulatory features to be displayed, such as Egr1. Be careful for the spelling and capitalisation of features. There are 3 cases. First, if you want to visualise only one feature (e.g. featureDisplay <- "CTCF") , you need to give the name of your specific feature. Second, if you want to visualise a set of features, you need to give a vector of features (e.g. featureDisplay <- c("Egr1","CTCF")) . Finally, if you want to visualise all features in the genomic region, use the word word "all" (e.g. featureDisplay <- "all"), "all" is set by default. You can find the complete list of features and their colour associated in annexe of vignette. Default value= "all"
}
  \item{datasetEnsembl}{
Allows the user to manually set which data set is used if required. Default value= "Null"
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
gen <- "hg19"
chr <- "chr7"
start <- 8000
end <- 80000
featureDisplay <- "CTCF"

if(interactive()){
 bindMotifsENSEMBLtrack<-BindingMotifsBiomart(gen,chr,start,end,featureDisplay)
  plotTracks(bindMotifsENSEMBLtrack, from = start, to = end) 
} else {
  data(bindMotifsENSEMBLtrack)
  plotTracks(bindMotifsENSEMBLtrack, from = start, to = end)
}

------------------------------------------------------------------------------

library("Gviz")
gen <- "hg19"
chr <- "chr7"
start <- 8000
end <- 80000
featureDisplay <- c("CTCF","Egr1")

if(interactive()){
 bindMotifsENSEMBLtrack<-BindingMotifsBiomart(gen,chr,start,end,featureDisplay)
  plotTracks(bindMotifsENSEMBLtrack, from = start, to = end) 
} else {
  data(bindMotifsENSEMBLtrack)
  plotTracks(bindMotifsENSEMBLtrack, from = start, to = end)
}

------------------------------------------------------------------------------
library("Gviz")
gen <- "hg19"
chr <- "chr7"
start <- 8000
end <- 80000
featureDisplay <- "all"
if(interactive()){
 bindMotifsENSEMBLtrack<-BindingMotifsBiomart(gen,chr,start,end,featureDisplay)
  plotTracks(bindMotifsENSEMBLtrack, from = start, to = end) 
} else {
  data(bindMotifsENSEMBLtrack)
  plotTracks(bindMotifsENSEMBLtrack, from = start, to = end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
