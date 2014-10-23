\name{genesNameENSEMBL}
\alias{genesNameENSEMBL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain the genes names in the genomic regions of interest from ENSEMBL
}
\description{
Obtain the genes names in the genomic regions of interest from ENSEMBL
}
\usage{
genesNameENSEMBL(gen, chr, start, end, dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{gen}{
the name of the genome
}
  \item{chr}{
the chromosome of interest
}
  \item{start}{
the first position in the region of interest (the smallest value)
}
  \item{end}{
the last position in the region of interest (the largest value)
}
  \item{dataset}{
 Name of the database to select genes
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
library("Gviz")
gen <- "hg19"
chr <- "chr7"
start <- 38290160
end <- 38303219

if(interactive()){
  dataset<- "hsapiens_gene_ensembl"
  geneNameEnsembl<- genesNameENSEMBL(gen,chr,start,end,dataset)
  geneNameEnsembl
} else {
  data(geneNameEnsembl)
  geneNameEnsembl
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

