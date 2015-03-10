\name{COSMICTrack}
\alias{COSMICTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one track of the genomic positions of variants from COSMIC
}
\description{
Create one track of the genomic positions of variants from COSMIC, the "Catalogue Of Somatic Mutations In Cancer" using the Gviz bioconductor package
}
\usage{
COSMICTrack(gen, chr, start, end, showId=FALSE)
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
http://bioconductor.org/packages/release/bioc/html/Gviz.html

http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=cosmic
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
\code{\link{CoreillCNVTrack}}, \code{\link{ClinVarMainTrack}}, 
\code{\link{ClinVarCnvTrack}}, 
}
\examples{
library("Gviz")
chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"
if(interactive()){
    cosmicVariant<-COSMICTrack(gen,chrom,start,end)
    plotTracks(cosmicVariant, from = start, to =end)
}else {
    data(cosmicVarianttrack)
    plotTracks(cosmicVariant, from = start, to =end)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
