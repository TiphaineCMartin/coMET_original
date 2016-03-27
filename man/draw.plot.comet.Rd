\name{draw.plot.comet}
\alias{draw.plot.comet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Display the three plots of coMET
}
\description{
Display the three plots of coMET according to configuration files : on
the upper plot is the plot of pvalue, the middle plot has the annotation
tracks, and the lower plot is the heatmap of correlation between CpG sites
}
\usage{
draw.plot.comet(config.var, gbl.var, newpage = TRUE)
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
\item{newpage}{
Option to ask if the plot should create on new page or not
}
}

\value{
Return gbl.var updated with the creation of elements composing the plot produced by the function comet
}

\author{
Tiphaine Martin
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ print }% __ONLY ONE__ keyword per line
