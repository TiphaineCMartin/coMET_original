\name{DNaseI_FANTOM}
\alias{DNaseI_FANTOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a enhancer/promoter track from FANTOM
}
\description{
Creates a track of promoters/enhancers from FANTOM using the Gviz bioconductor package. A complete list of features and their associated colours can be found in the user guide.
}
\usage{
DNaseI_FANTOM(gen="hg19", chr, start, end, bedFilePath, featureDisplay='enhancer', stacking_type="dense")
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
A vector of regulatory features to be displayed, such as enhancer. Spelling and capitalisation of features must be identical to those in the user guide. There are three possibilities. First, the visualisation of only one feature (e.g. featureDisplay <- "Predicted heterochomatin"), only the name of the specific feature is required. Second, visualisation of a set of features, for this a vector of features is required (e.g. featureDisplay <- c("enhancer","promoter")). Finally, visualison all features in the genomic region, achived by using the word "all" (e.g. featureDisplay <- "all"), "all" is set by default. You can find the complete list of features and their associated colours in the user guide.
}
  \item{stacking_type}{
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
chr<- "chr1"
start <- 6000000
end <- 6500000

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
enhFantomFile <- file.path(extdata, "/FANTOM/human_permissive_enhancers_phase_1_and_2.bed")

if(interactive()){
 enhFANTOMtrack <- DNaseI_FANTOM(gen,chr,start, end, enhFantomFile, featureDisplay='enhancer')
  plotTracks(enhFANTOMtrack, from = start, to = end)
} else {
  data(enhFANTOMtrack)
  plotTracks(enhFANTOMtrack, from = start, to = end)
}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
