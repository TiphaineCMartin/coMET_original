\name{TFBS_FANTOM}
\alias{TFBS_FANTOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a TFBS motif track from FANTOM
}
\description{
Creates a track of TFBS motifs from FANTOM using the Gviz bioconductor package. A complete list of features and their associated colours can be found in the user guide.
}
\usage{
TFBS_FANTOM(gen, chr, start, end, bedFilePath)
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
The path of the BED file from Kheradpour and Kellis, 2014.
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
AP1FantomFile <- file.path(extdata, "/FANTOM/Fantom_hg19.AP1_MA0099.2.sites.txt")

if(interactive()){
  tfbsFANTOMtrack <- TFBS_FANTOM(gen,chr,start, end, AP1FantomFile)
  plotTracks(tfbsFANTOMtrack, from = start, to = end)
} else {
  data(tfbsFANTOMtrack)
  plotTracks(tfbsFANTOMtrack, from = start, to = end)
}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
