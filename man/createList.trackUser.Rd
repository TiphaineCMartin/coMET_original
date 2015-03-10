\name{createList.trackUser}
\alias{createList.trackUser}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create list of Gviz's tracks from user's data
}
\description{
Create list of Gviz's tracks from user's data
}
\usage{
createList.trackUser(config.var, gbl.var)
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
Update the object listtracks_user of gbl.var from the user data if the user gives some data for the annotation tracks.
}

\author{
Tiphaine Martin
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ dplot }
