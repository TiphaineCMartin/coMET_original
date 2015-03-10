\name{chrUCSC2ENSEMBL}
\alias{chrUCSC2ENSEMBL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Removing "chr" to the chromosome number from UCSC to transform it to ENSEMBL chromosome format
}
\description{
Removing "chr" at the beginning of the chromosome number
}
\usage{
chrUCSC2ENSEMBL(chr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{chr}{
the chromosome number in UCSC format 
}
}

\value{
    the number of chromosome at ENSEMBL format
}

\author{
Tiphaine Martin
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
    chr<-"chr7"
    chrUCSC2ENSEMBL(chr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dplot }
