\name{snpBiomart}
\alias{snpBiomart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a short variation track from ENSEMBL
}
\description{
Create a 'Short Variation' track from ENSEMBL using the Gviz bioconductor package
}
\usage{
snpBiomart(chr, start, end, dataset, showId=FALSE, title = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
  \item{title}{
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
 \code{\link{snpLocationsUCSC}},  \code{\link{structureBiomart}},  \code{\link{COSMICTrack}}, 
  \code{\link{CoreillCNVTrack}}, \code{\link{ClinVarMainTrack}}, 
  \code{\link{ClinVarCnvTrack}}, 
}
\examples{
library("Gviz")
gen <- "hg19"
chr <- "chr2"
start <- 38290160
end <- 38303219

if(interactive()){
  snptrack <- snpBiomart(chr, start, end, 
                       dataset="hsapiens_snp_som",showId=FALSE)
  plotTracks(snptrack, from=start, to=end) 
} else {
  data(snpBiomarttrack)
  plotTracks(snptrack, from=start, to=end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
