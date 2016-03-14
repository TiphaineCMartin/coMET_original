\name{interestTranscript_ENSEMBL}
\alias{interestTranscript_ENSEMBL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a track of transcripts from ENSEMBL
}
\description{
Create a track to visualize different transcripts from ENSEMBL using the Gviz bioconductor package
}
\usage{
interestTranscript_ENSEMBL(gen, chr, start, end,interestfeatures,interestcolor,showId = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gen}{
the name of the genome
}
  \item{chr}{
the chromosome of interest
}
  \item{start}{
the first position in the region of interest (the smallest value)
}
  \item{end}{
the last position in the region of interest (the largest value)
}
  \item{interestfeatures}{
A data frame with 3 columns: start of features, end of features, and type of features
}
  \item{interestcolor}{
A list with the color for each new features defined
}
  \item{showId}{
  Show the ID of the genetic elements
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
A BiomartGeneRegionTrack object of Gviz
}
\references{
http://bioconductor.org/packages/release/bioc/html/Gviz.html

http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=ensGene
}
\author{
Tiphaine Martin
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{ISCA_UCSC}},  \code{\link{GWAScatalog_UCSC}}, \code{\link{knownGenes_UCSC}},
 \code{\link{genesName_ENSEMBL}},  \code{\link{GeneReviews_UCSC}}, \code{\link{GAD_UCSC}},  \code{\link{genes_ENSEMBL}}, \code{\link{xenorefGenes_UCSC}},
}
\examples{
library("Gviz")
gen <- "hg19"
chr <- "chr15"
start <- 75011669
end <- 75019876
interestfeatures <- rbind(c("75017782","75017835","bad"),c("75013755","75013844","good"))
interestcolor <- list("bad"="red", "good"="green")

if(interactive()){
  interesttransENSMBLtrack<-interestTranscript_ENSEMBL(gen,chr,start,end,interestfeatures,interestcolor,showId=TRUE)
  plotTracks(interesttransENSMBLtrack, from=start, to=end)
} else {
  data(interesttransENSMBLtrack)
  plotTracks(interesttransENSMBLtrack, from=start, to=end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
