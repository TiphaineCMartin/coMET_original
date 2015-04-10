\name{coMET-package}
\alias{coMET-package}
\alias{coMET}
\docType{package}
\title{
visualisation of regional epigenome-wide association scan (EWAS) results and DNA co-methylation patterns (and also for other omic-WAS)
}
\description{
coMET is an R package for visualising EWAS results in a genomic region. Along with phenotype-association plots, coMET also generates plots of co-methylation patterns and provides a series of annotation tracks. The software is designed for epigenetic data, but can also be applied to genomic and functional genomic datasets (other omic-WAS results) in any species.
}
\details{
\tabular{ll}{
Package: \tab coMET\cr
Type: \tab Package\cr
Version: \tab 0.99.10\cr
Date: \tab 2015-04-10\cr
License: \tab GPL (>=2)\cr
}
coMET is an R package that can generate regional plots of EWAS results, DNA co-methylation patterns, and genomic information. A coMET figure includes 3 panels with a plot of P-values from EWAS, customized annotation tracks, and a triangle heatmap plot which demonstrates the correlation structure of DNA methylation at the CpG sites in the genomic region. Plots are created as PDF or EPS files.
}
\author{
Tiphaine C. Martin, Idil Yet, Pei-Chien Tsai, Jordana T. Bell

Maintainer: Tiphaine Martin <tiphaine.martin@kcl.ac.uk>

Website: http://www.epigen.kcl.ac.uk/comet
}
\references{
Martin, T.C, Erte, I, Tsai, P-C, Bell, J.T., coMET: an R plotting package to visualize regional plots of epigenome-wide association scan results, QG14, 2014.
}
\keyword{ package }
%\seealso{
% \code{\link[<coMET>:<coMET>-package]{<coMET>}}
%}
\examples{

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4comet.txt") 
myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")

chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"

if(interactive()){
    genetrack <-genesENSEMBL(gen,chrom,start,end,showId=TRUE)
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
                    
    comet(config.file=configfile, mydata.file=myinfofile,  mydata.type="listfile",
         cormatrix.file=mycorrelation, cormatrix.type="listfile",
         mydata.file=myexpressfile, mydata.large.type="listfile",
         tracks.gviz=listgviz, 
         verbose=FALSE, print.image=FALSE,disp.pvalueplot=TRUE)
} else {
    data(geneENSEMBLtrack)
    data(snpBiomarttrack)
    data(ISCAtrack)
    data(strucBiomarttrack)
    data(ClinVarCnvTrack)
    data(clinVarMaintrack)
    data(GWASTrack)
    data(GeneReviewTrack)

    listgviz <- list(genetrack,snptrack,strutrack,clinVariant,
                    clinCNV,gwastrack,geneRtrack)
    comet(config.file=configfile, mydata.file=myinfofile, mydata.type="listfile",
        cormatrix.file=mycorrelation, cormatrix.type="listfile",
        mydata.large.file=myexpressfile, mydata.large.type="listfile",
        tracks.gviz=listgviz, 
        verbose=FALSE, print.image=FALSE,disp.pvalueplot=TRUE)
}
}
