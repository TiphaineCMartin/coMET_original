\name{HistoneOne}
\alias{HistoneOne}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create track one type of Histone density from UCSC
}
\description{
Create track one type of Histone density from UCSC (ENCODE/Broad)
}
\usage{
HistoneOne(gen, chr, start, end, mySession, track.name = NULL, table.name = NULL)
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
  \item{mySession}{
the object session from the function browserSession of rtracklayer
}
 \item{track.name}{
the name of track, for example : Histone
}
  \item{table.name}{
the name of table from the track
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
