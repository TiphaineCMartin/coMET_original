\name{CoreillCNVTrack}
\alias{CoreillCNVTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the genomic positions of copy-number variants (CNVs) in chromosomal aberration and inherited disorder cell lines in the NIGMS Human Genetic Cell Repository
}
\description{
Get the genomic positions of copy-number variants (CNVs) in chromosomal aberration and inherited disorder cell lines in the NIGMS Human Genetic Cell Repository
}
\usage{
CoreillCNVTrack(gen, chr, start, end, showId=FALSE)
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
 \item{showId}{
Show Id of genetic elements
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

coreilVariant<-CoreillCNVTrack(gen,chrom,start,end)
plotTracks(coreilVariant)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
