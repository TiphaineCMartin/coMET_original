\name{structureBiomart}
\alias{structureBiomart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create track Structural Variation from ENSEMBL
}
\description{
Create track Structural Variation from ENSEMBL
}
\usage{
structureBiomart(chr, start, end, strand, dataset, showId=FALSE, title = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{
the chromosome of interest
}
  \item{start}{
the first position of region of interest (the smallest value)
}
  \item{end}{
the last position of region of interest (the biggest value)
}
  
  \item{strand}{
 the strand of genome to extract structure
}
  \item{dataset}{
 The name of database
}
  \item{showId}{
  Show the id of element or not
}
  \item{title}{
 The name of annotation track
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
