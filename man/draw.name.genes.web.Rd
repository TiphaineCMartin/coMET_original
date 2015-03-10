\name{draw.name.genes.web}
\alias{draw.name.genes.web}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
display the gene names 
}
\description{
display the gene names for the web page(see cometweb)
}
\usage{
draw.name.genes.web(config.var, gbl.var)
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
Updated gbl.var with the list of names of genes found the region of interest. This function is called only in comet.web
}

\author{
Tiphaine Martin
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ print }% __ONLY ONE__ keyword per line
