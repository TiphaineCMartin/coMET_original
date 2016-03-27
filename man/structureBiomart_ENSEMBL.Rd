\name{structureBiomart_ENSEMBL}
\alias{structureBiomart_ENSEMBL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a structural variation track from ENSEMBL
}
\description{
Create a 'Structural Variation' track from ENSEMBL using the Gviz bioconductor package
}
\usage{
structureBiomart_ENSEMBL(gen, chr, start, end, strand, dataset, showId=FALSE, title_track = NULL)
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

  \item{strand}{
 the strand to extract structure data for
}
  \item{dataset}{
 The name of the database. Example "hsapiens_structvar_som"
}
  \item{showId}{
  Show the the ID of the element
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
 \code{\link{snpLocations_UCSC}},  \code{\link{snpBiomart_ENSEMBL}},  \code{\link{COSMIC_UCSC}},
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
  strutrack <- structureBiomart_ENSEMBL(chr, start, end,
                              strand, dataset="hsapiens_structvar_som")
  plotTracks(strutrack, from=start, to=end)
}else {
  data(strucBiomarttrack)
  plotTracks(strutrack, from=start, to=end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
