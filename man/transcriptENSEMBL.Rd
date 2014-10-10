\name{transcriptENSEMBL}
\alias{transcriptENSEMBL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creation track for transcription sites
}
\description{
  Creation track for transcription sites in using the Gviz bioconductor package
}
\usage{
transcriptENSEMBL(gen, chr, start, end,showId = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gen}{
the name of genome
}
  \item{chr}{
the chromosome of interest
}
  \item{start}{
the first position of region of interest (the smallest value)
}
  \item{end}{
the last position of region of interest (the biggest value)
}
  \item{showId}{
  Show Id of genetic elements
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A BiomartGeneRegionTrack object of Gviz
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tiphaine Martin
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library("Gviz")
gen <- "hg19"
chr <- "chr7"
start <- 38290160
end <- 38303219
tmp<-transcriptENSEMBL(gen,chr,start,end)
plotTracks(tmp, from = start, to =end) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
