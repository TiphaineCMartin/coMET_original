\name{ClinVarCnvTrack}
\alias{ClinVarCnvTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one track of the genomic positions of variants from the ClinVar database (CNV only)
}
\description{
Create one track of the genomic positions of variants from the ClinVar database (CNV only, Variants excluded) using the Gviz bioconductor package
}
\usage{
ClinVarCnvTrack(gen, chr, start, end, showId = FALSE)
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
    \code{\link{snpLocationsUCSC}}, \code{\link{structureBiomart}}, 
    \code{\link{snpBiomart}}, 
    \code{\link{CoreillCNVTrack}}, \code{\link{COSMICTrack}}, 
    \code{\link{ClinVarMainTrack}}
}
\examples{
library("Gviz")
chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"
if(interactive()){
    clinCNV<-ClinVarCnvTrack(gen,chrom,start,end)
    plotTracks(clinCNV, from = start, to =end) 
}else {
    data(ClinVarCnvTrack)
    plotTracks(clinCNV, from = start, to =end) 
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

