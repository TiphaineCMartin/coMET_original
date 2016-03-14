\name{create.color.bar}
\alias{create.color.bar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create color bar of heatmap
}
\description{
Create color legend for the correlation matrix
}
\usage{
create.color.bar(config.var, gbl.var)
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
}

\value{
 list of different matrix
  \item{color.cut}{the matrix of color related to correlation matrix} ,
  \item{color.cut.ref}{the vector of color related to the reference CpG} ,
  \item{cormatrix.key}{the generic panel having different panels associated with the correlation matrix},
  \item{map.label.ldtype}{the panel with the method of creation of correlation matrix}.
  \item{map.label.distance}{the panel with the legend of distance}
}

\author{
Tiphaine Martin
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
