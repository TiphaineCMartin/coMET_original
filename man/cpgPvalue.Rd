\name{cpgPvalue}
\alias{cpgPvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a plot of pvalue of CpG with DataTrack of Gviz
}
\description{
Create a plot of pvalue of CpG with DataTrack of Gviz
}
\usage{
cpgPvalue(cprange, data, chr, start, end, typefunction, title)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{cprange}{
Range to visualise
}
    \item{data}{
data to analyse
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
    \item{typefunction}{
    Type of function to visualise the data
}
    \item{title}{
    Name of tracks
}
}

\value{
the object DataTrack of Gviz
}

\author{
Tiphaine Martin
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

