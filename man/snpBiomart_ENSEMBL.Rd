\name{snpBiomart_ENSEMBL}
\alias{snpBiomart_ENSEMBL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a short variation track from ENSEMBL
}
\description{
Create a 'Short Variation' track from ENSEMBL using the Gviz bioconductor package
}
\usage{
snpBiomart_ENSEMBL(gen,chr, start, end, dataset, showId=FALSE, title_track = NULL)
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
  \item{dataset}{
  The name of the database. Example "hsapiens_snp_som"
}
  \item{showId}{
  Show the the ID of element or not
}
  \item{title_track}{
 The name of the annotation track
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
An AnnotationTrack object of Gviz
}
\references{
Go to ENSEMBL Biomart

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
chr <- "chr2"
start <- 38290160
end <- 38303219

if(interactive()){
  snptrack <- snpBiomart_ENSEMBL(gen,chr, start, end,
                       dataset="hsapiens_snp",showId=FALSE)
  plotTracks(snptrack, from=start, to=end)
} else {
  data(snpBiomarttrack)
  plotTracks(snptrack, from=start, to=end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
