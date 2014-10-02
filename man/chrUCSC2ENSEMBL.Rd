\name{chrUCSC2ENSEMBL}
\alias{chrUCSC2ENSEMBL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adding "chr" from the chromsome of UCSC to become ENSEMBL's chromosome
}
\description{
Adding the letter chr in the beginning of the name of chromosome from UCSC
}
\usage{
chrUCSC2ENSEMBL(chr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{chr}{
the name of chromosome at UCSC format 
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
chr<-"chr7"
chrUCSC2ENSEMBL(chr)


## The function is currently defined as
function (chr) 
{
    gsub("^chr", "", chr)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
