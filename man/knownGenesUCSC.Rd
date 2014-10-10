\name{knownGenesUCSC}
\alias{knownGenesUCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create track Known genes from UCSC
}
\description{
create track Known genes from UCSC in using the Gviz bioconductor package
}
\usage{
knownGenesUCSC(gen, chr, start, end)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gen}{
the name of genome
}
  \item{chr}{
the chromosome of interest
}
  \item{start}{
the first position of region of interest (the smallest value)
}
  \item{end}{
the last position of region of interest (the biggest value)
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
chr <- "chr7"
start <- 38290160
end <- 38303219
tmp<-knownGenesUCSC(gen,chr,start,end)
plotTracks(tmp, from = start, to =end) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
