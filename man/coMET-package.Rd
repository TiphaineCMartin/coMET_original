\name{coMET-package}
\alias{coMET-package}
\alias{coMET}
\docType{package}
\title{
coMET: Visualisation of pvalue and correlation between genomic data such as DNA methylation
}
\description{
 The coMET is a R package to visualize the EWAS (epigenome-wide association scans) results in a genomic region. coMET package generates the plots of association, comethylation patterns and a series of annotation tracks at genomic scale.
}
\details{
\tabular{ll}{
Package: \tab coMET\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2014-09-26\cr
License: \tab GPL (>=2)\cr
}
coMET package that can generate the regional plot capturing the features of co-methylation patterns, EWAS results, and genomic information. A coMET figure includes plot of p-value from the EWAS result, provides customized annotation tracks, and a triangle heatmap plot which demonstrates the correlation structure of CpG sites in a genomic region, calculated by the pairwise Spearman’s rank correlation method.Plots are created as PDF,EPS files.

A list containing two items: config.var and gbl.var, which includes the values of all significant variables used by coMET.
}
\author{
Tiphaine Martin, Idil Erte, Pei-Chien Tsai, Jordana T. Bell

Maintainer: Tiphaine Martin <tiphaine.martin@kcl.ac.uk>
Website: http://www.epigen.kcl.ac.uk/comet
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
~~ documentation directory ~~
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<coMET>:<coMET>-package]{<coMET>}} ~~
}
\examples{
library(coMET)
library(Gviz)
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4comet.txt") 
chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"

genetrack <-genesENSEMBL(gen,chrom,start,end,showId=FALSE)
snptrack <- snpBiomart(chrom, start, end, 
                       dataset="hsapiens_snp_som",showId=FALSE)
strutrack <- structureBiomart(chrom, start, end, 
                              strand, dataset="hsapiens_structvar_som")
clinVariant<-ClinVarMainTrack(gen,chrom,start,end)
clinCNV<-ClinVarCnvTrack(gen,chrom,start,end)
gwastrack <-GWASTrack(gen,chrom,start,end)
geneRtrack <-GeneReviewsTrack(gen,chrom,start,end)

listgviz <- list(genetrack,snptrack,strutrack,clinVariant,
                 clinCNV,gwastrack,geneRtrack)
comet(config.file=configfile,TRACKS.GVIZ=listgviz, 
      VERBOSE=FALSE, PRINT.IMAGE=FALSE,DISP.PVALUEPLOT=FALSE)
}
