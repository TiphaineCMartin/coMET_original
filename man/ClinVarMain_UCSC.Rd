\name{ClinVarMain_UCSC}
\alias{ClinVarMain_UCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one track of the genomic positions of variants from the ClinVar database (variants only)
}
\description{
Create one track of the genomic positions of variants from the ClinVar database (Variants only, CNV excluded) using the Gviz bioconductor package
}
\usage{
ClinVarMain_UCSC(gen, chr, start, end, showId=FALSE)
}

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
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
An UcscTrack object of Gviz
}
\references{
http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=clinvar

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
    \code{\link{snpLocations_UCSC}}, \code{\link{structureBiomart_ENSEMBL}},
    \code{\link{snpBiomart_ENSEMBL}},
    \code{\link{CoreillCNV_UCSC}}, \code{\link{COSMIC_UCSC}},
    \code{\link{ClinVarCnv_UCSC}}
}
\examples{
library("Gviz")
gen <- "hg38"
chrom <- "chr2"
start <- 100000
end <- 10000000

if(interactive()) {
    clinVariant<-ClinVarMain_UCSC(gen,chrom,start,end)
    plotTracks(clinVariant, from = start, to =end)
}else{
    data(clinVarMaintrack)
    plotTracks(clinVariant, from = start, to =end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
