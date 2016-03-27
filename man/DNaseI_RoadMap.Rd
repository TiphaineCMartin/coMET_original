\name{DNaseI_RoadMap}
\alias{DNaseI_RoadMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a promoter/enhancer regions track from a file of RoadMap
}
\description{
Creates a track of promoter/enhancer regions from a file of RoadMap using the Gviz bioconductor package. A complete list of features and their associated colours can be found in the user guide.
}
\usage{
DNaseI_RoadMap(gen="hg19", chr, start, end, bedFilePath, featureDisplay='promotor',showId=TRUE, type_stacking="dense")
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
  \item{showId}{
Allows to visualise the Id of DNAse group.
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
start <- 707612
end <- 722151
gen="hg19"

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
bedFilePath <- file.path(extdata, "RoadMap/regions_prom_E063.bed")

if(interactive()){
  DNaseI_RoadMapSingle <- DNaseI_RoadMap(gen,chr,start, end, bedFilePath, featureDisplay='promotor'  )
  plotTracks(DNaseI_RoadMapSingle, from = start, to = end)
} else {
  data(DNaseI_RoadMapSingle)
  plotTracks(DNaseI_RoadMapSingle, from = start, to = end)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
