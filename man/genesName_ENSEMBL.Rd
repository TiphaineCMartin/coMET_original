\name{genesName_ENSEMBL}
\alias{genesName_ENSEMBL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain the genes names in the genomic regions of interest from ENSEMBL
}
\description{
Obtain the genes names in the genomic regions of interest from ENSEMBL
}
\usage{
genesName_ENSEMBL(gen, chr, start, end, dataset)
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
Can be null
}
\value{
List of name of genes found in this region of interest.
}
\references{
go to ENSEMBL

http://bioconductor.org/packages/release/bioc/html/Gviz.html
}
\author{
Tiphaine Martin
}

\seealso{
 \code{\link{ISCA_UCSC}},  \code{\link{GWAScatalog_UCSC}}, \code{\link{knownGenes_UCSC}},
 \code{\link{GeneReviews_UCSC}},  \code{\link{GAD_UCSC}}, \code{\link{genes_ENSEMBL}},  \code{\link{xenorefGenes_UCSC}}, \code{\link{transcript_ENSEMBL}},
}
\examples{
library("Gviz")
gen <- "hg38"
chr <- "chr7"
start <- 38290160
end <- 38303219

if(interactive()){
  dataset<- "hsapiens_gene_ensembl"
  geneNameEnsembl<- genesName_ENSEMBL(gen,chr,start,end,dataset)
  geneNameEnsembl
} else {
  data(geneNameEnsembl)
  geneNameEnsembl
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

