\name{cpgIslandsUCSC}
\alias{cpgIslandsUCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the correlation matrix between CpG sites
}
\description{
Compute the correlation matrix between CpG sites
}
\usage{
cpgIslandsUCSC(gen, chr, start, end)
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
chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"

tmp<-cpgIslandsUCSC(gen, chr, start, end)
plotTracks(tmp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~method }
\keyword{ ~hplot }% __ONLY ONE__ keyword per line
