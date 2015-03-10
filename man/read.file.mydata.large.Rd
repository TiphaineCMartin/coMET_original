\name{read.file.mydata.large}
\alias{read.file.mydata.large}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read the files of supplemantary data and extract data
}
\description{
Read the files of supplemantary data and extract data
}
\usage{
read.file.mydata.large(large.split.mydata.file, config.var, gbl.var, numfile.large)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{large.split.mydata.file}{
List of supplementary files to check format
}
\item{config.var}{
list of all variables defined in configuration file or via options of
comet function
}
\item{gbl.var}{
list of internal variables
}
\item{numfile.large}{
The number of file to check in the list
}
}
\value{
Return gbl.var updated after reading the extra data 
}
\author{
Tiphaine Martin
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ dplot }% __ONLY ONE__ keyword per line
