\name{ClinVarMainTrack}
\alias{ClinVarMainTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one track of the genomic positions of variants from the ClinVar database (variants only)
}
\description{
Create one track of the genomic positions of variants from the ClinVar database (Variants only, CNV excluded) using the Gviz bioconductor package
}
\usage{
ClinVarMainTrack(gen, chr, start, end, showId=FALSE)
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
the first position in region of interest (the smallest value)
}
  \item{end}{
the last position in region of interest (the biggest value)
}
 \item{showId}{
Show the ID of the genetic elements
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An UcscTrack object of Gviz
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
chrom <- "chr2"
start <- 38290160
end <- 38303219

if(interactive()) {
  clinVariant<-ClinVarMainTrack(gen,chrom,start,end)
  plotTracks(clinVariant, from = start, to =end) 
}else{
  data(clinVarMaintrack)
  plotTracks(clinVariant, from = start, to =end) 
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
