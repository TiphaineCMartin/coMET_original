\name{cpgIslandsUCSC}
\alias{cpgIslandsUCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create track CpG Island from UCSC
}
\description{
create track CpG Island from UCSC using the Gviz bioconductor package
}
\usage{
cpgIslandsUCSC(gen, chr, start, end)
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
chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"

if(interactive()) {
  cpgIstrack<-cpgIslandsUCSC(gen, chrom, start, end)
  plotTracks(cpgIstrack, from = start, to =end)
}else {
  data(cpgIslandtrack)
  plotTracks(cpgIstrack, from = start, to =end)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
