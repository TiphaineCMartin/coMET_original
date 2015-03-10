\name{read.file.cormatrix}
\alias{read.file.cormatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read, compute and extract the values from correlation matrix file
}
\description{
Read, compute and extract the values from correlation matrix file
}
\usage{
read.file.cormatrix(config.var, gbl.var, split.cormatrix.file = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{config.var}{
list of all variables defined in configuration file or via options of
comet function
}
\item{gbl.var}{
list of internal variables
}
\item{split.cormatrix.file}{
File of correlation matrix 
}
}
\value{
Return gbl.var updated with the raw data of the correlation matrix and the correlation matrix computed
}
\author{
Tiphaine Martin
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ dplot }% __ONLY ONE__ keyword per line
