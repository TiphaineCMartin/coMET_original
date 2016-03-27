\name{refGenes_UCSC}
\alias{refGenes_UCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a track of RefSeq genes from the UCSC genome browser
}
\description{
Create a track of RefSeq genes from the UCSC genome browser using the Gviz bioconductor package
}
\usage{
refGenes_UCSC(gen, chr, start, end, IdType="Ref", showId=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gen}{
The name of the genome
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
  \item{IdType}{
When set to 'ref' shows the gene reference, when set to "name" shows the gene name
}
 \item{showId}{
Shows the ID or name of the genetic elements
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
An UcscTrack object of Gviz
}
\references{
http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=knownGene

http://bioconductor.org/packages/release/bioc/html/Gviz.html
}
\author{
Tiphaine Martin

Tom Hardiman
}
%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{ISCA_UCSC}}, \code{\link{GWAScatalog_UCSC}}, \code{\link{genesName_ENSEMBL}},
 \code{\link{GeneReviews_UCSC}},  \code{\link{GAD_UCSC}}, \code{\link{genes_ENSEMBL}},  \code{\link{xenorefGenes_UCSC}}, \code{\link{transcript_ENSEMBL}}, \code{\link{knownGenes_UCSC}}
}
\examples{
library("Gviz")
gen <- "hg38"
chr <- "chr7"
start <- 38203219
end <- 38303219
IdType <- "name"

if(interactive()) {
  genesUcsctrack<-refGenes_UCSC(gen,chr,start,end,IdType)
  plotTracks(genesUcsctrack, from = start, to =end)
}else {
  data(genesUcsctrack)
  plotTracks(genesUcsctrack, from = start, to =end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
