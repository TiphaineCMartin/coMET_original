\name{ISCATrack}
\alias{ISCATrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one track of the genomic positions of variants from ISCA
}
\description{
Create one track of the genomic positions of variants from International Standards for Cytogenomic Arrays (ISCA) Consortium using the Gviz bioconductor package
}
\usage{
ISCATrack(gen, chr, start, end, mySession, table.name, showId=FALSE)
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
 \item{mySession}{
the object session from the function browserSession of rtracklayer
}
  \item{table.name}{
A table of ISCAT classifications: iscaBenign, iscaCuratedBenign, iscaCuratedPathogenic, iscaLikelyBenign, iscaLikelyPathogenic, iscaPathGainCum, iscaPathLossCum, iscaPathogenic, iscaUncertain
}
 \item{showId}{
Show the ID of the genetic elements
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
An UcscTrack object of Gviz
}
\references{
http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=iscaComposite

http://bioconductor.org/packages/release/bioc/html/Gviz.html
}
\author{
Tiphaine Martin
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{GWASTrack}},  \code{\link{knownGenesUCSC}}, \code{\link{genesNameENSEMBL}}, 
 \code{\link{GeneReviewsTrack}},  \code{\link{GADTrack}}, \code{\link{genesENSEMBL}},  \code{\link{xenorefGenesUCSC}}, \code{\link{transcriptENSEMBL}},
}
\examples{
library("Gviz")
gen <- "hg19"
chr <- "chr2"
start <- 38292433
end <- 38305492

if(interactive()){
  BROWSER.SESSION="UCSC"
  mySession <- browserSession(BROWSER.SESSION)
  genome(mySession) <- gen
  iscatrack <-ISCATrack(gen,chrom,start,end,mySession, table="iscaPathogenic")
  plotTracks(iscatrack, from = start, to =end)
} else {
  data(ISCAtrack)
  plotTracks(iscatrack, from = start, to =end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
