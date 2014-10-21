### R code from vignette source 'vignettes/coMET/inst/doc/coMET.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: style knitr
###################################################
library(knitr)
opts_chunk$set(fig.path='figure/minimal-', fig.align='center', fig.show='hold')
options(replace.assign=TRUE,width=90)


###################################################
### code chunk number 3: install coMET
###################################################
source("http://bioconductor.org/biocLite.R")
biocLite("coMET")


###################################################
### code chunk number 4: load coMET
###################################################
library(coMET)

###################################################
### code chunk number 5: help
###################################################
?comet
?comet.web

###################################################
### code chunk number 6: read info file (site without association)
###################################################
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
infofile <- file.path(extdata, "cyp1b1_infofile.txt")

data_info <-read.csv(infofile, header = TRUE,
                     sep = "\t", quote = "")

head(data_info)

###################################################
### code chunk number 7: read supplementary file (region with association)
###################################################
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
infoexp <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
#infoexp <- "../inst/extdata/cyp1b1_infofile_exprGene_region.txt" 

data_infoexp <-read.csv(infoexp, header = TRUE,
                        sep = "\t", quote = "")

head(data_infoexp)

###################################################
### code chunk number 8: read correlation matrix
###################################################
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
corfile <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
#corfile <- "../inst/extdata/cyp1b1_res37_rawMatrix.txt" 

data_cor <-read.csv(corfile, header = TRUE,
                    sep = "\t", quote = "")
data_cor[1:6,1:6]


###################################################
### code chunk number 9: read configuration file
###################################################
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4webserver.txt")
#configfile <- "../inst/extdata/config_cyp1b1_zoom_4webserver.txt" 

data_config <-read.csv(configfile, quote = "")
data_config

###################################################
### code chunk number 10: Creation plot with comet.web
###################################################
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt") 
configfile <- file.path(extdata, "config_cyp1b1_zoom_4webserver.txt")
comet.web(config.file=configfile, MYDATA.FILE=myinfofile, 
          CORMATRIX.FILE=mycorrelation ,MYDATA.LARGE.FILE=myexpressfile, 
          PRINT.IMAGE=FALSE,VERBOSE=FALSE)


###################################################
### code chunk number 11: Creation plot with comet
###################################################
library(Gviz)
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4comet.txt")
myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
#configfile <- "../inst/extdata/config_cyp1b1_zoom_4comet.txt" 
chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"


if(interactive()) {
  BROWSER.SESSION="UCSC"
  mySession <- browserSession(BROWSER.SESSION)
  genome(mySession) <- gen
  
  genetrack <-genesENSEMBL(gen,chrom,start,end,showId=FALSE)
  snptrack <- snpBiomart(chrom, start, end, dataset="hsapiens_snp_som",showId=FALSE)
  strutrack <- structureBiomart(chrom, start, end, strand,
                                dataset="hsapiens_structvar_som",showId=FALSE)
  iscatrack <-ISCATrack(gen,chrom,start,end,mySession, table="iscaPathogenic")
  
  listgviz <- list(genetrack,snptrack,iscatrack)
  comet(config.file=configfile, MYDATA.FILE=myinfofile, CORMATRIX.FILE=mycorrelation,
      MYDATA.LARGE.FILE=myexpressfile, TRACKS.GVIZ=listgviz, 
      VERBOSE=FALSE, PRINT.IMAGE=FALSE)
} else {
  data(geneENSEMBLtrack)
  data(snpBiomarttrack)
  data(ISCAtrack)
  
  listgviz <- list(genetrack,snptrack,iscatrack)
  comet(config.file=configfile, MYDATA.FILE=myinfofile, CORMATRIX.FILE=mycorrelation, 
      MYDATA.LARGE.FILE=myexpressfile, TRACKS.GVIZ=listgviz, 
      VERBOSE=FALSE, PRINT.IMAGE=FALSE)
}


###################################################
### code chunk number 11: Creation plot with comet without pvalue plot
###################################################
library(Gviz)
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4comet.txt")
#configfile <- "../inst/extdata/config_cyp1b1_zoom_4comet.txt" 
chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"

if(interactive()){
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
} else {
  data(geneENSEMBLtrack)
  data(snpBiomarttrack)
  data(strucBiomarttrack)
  data(ClinVarCnvTrack)
  data(clinVarMaintrack)
  data(GWASTrack)
  data(GeneReviewTrack)
  
  listgviz <- list(genetrack,snptrack,strutrack,clinVariant,
                   clinCNV,gwastrack,geneRtrack)
  comet(config.file=configfile,TRACKS.GVIZ=listgviz, 
        VERBOSE=FALSE, PRINT.IMAGE=FALSE,DISP.PVALUEPLOT=FALSE)
}



