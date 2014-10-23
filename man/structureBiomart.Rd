\name{structureBiomart}
\alias{structureBiomart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a structural variation track from ENSEMBL
}
\description{
Create a 'Structural Variation' track from ENSEMBL using the Gviz bioconductor package
}
\usage{
structureBiomart(chr, start, end, strand, dataset, showId=FALSE, title = NULL)
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
  
  \item{strand}{
 the strand to extract structure data for
}
  \item{dataset}{
 The name of the database. Example "hsapiens_structvar_som"
}
  \item{showId}{
  Show the the ID of the element
}
  \item{title}{
 The name of the annotation track
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An AnnotationTrack object of Gviz
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
chr <- "chr2"
start <- 38290160
end <- 38303219

if(interactive()){
  strutrack <- structureBiomart(chrom, start, end, 
                              strand, dataset="hsapiens_structvar_som")
  plotTracks(strutrack, from = start, to =end) 
}else {
  data(strucBiomarttrack)
  plotTracks(strutrack, from = start, to =end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
