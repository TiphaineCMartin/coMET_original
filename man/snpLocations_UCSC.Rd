\name{snpLocations_UCSC}
\alias{snpLocations_UCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a SNP track from UCSC
}
\description{
Create a SNP track from UCSC using the Gviz bioconductor package
}
\usage{
snpLocations_UCSC(gen, chr, start, end, track="All SNPs(142)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{gen}{
the name of the genome. Data is not currently available for GRCh38 (hg38).
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

  \item{track}{
 The name of the database. Default "All SNPs(142)"
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
An UcscTrack object of Gviz
}
\references{
http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=snp141

http://bioconductor.org/packages/release/bioc/html/Gviz.html
}
\author{
Tiphaine Martin
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{snpLocations_UCSC}},  \code{\link{structureBiomart_ENSEMBL}},  \code{\link{COSMIC_UCSC}},
  \code{\link{CoreillCNV_UCSC}}, \code{\link{ClinVarMain_UCSC}},
  \code{\link{ClinVarCnv_UCSC}},
}
\examples{
library("Gviz")
gen <- "hg38"
chr <- "chr7"
start <- 38290160
end <- 38303219

if(interactive()) {
  snpUCSCtrack<-snpLocations_UCSC(gen,chr,start,end,"All SNPs(142)")
  plotTracks(snpUCSCtrack, from = start, to =end)
} else {
  data(snpUCSCtrack)
  plotTracks(snpUCSCtrack, from = start, to =end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
