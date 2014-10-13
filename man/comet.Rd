\name{comet}
\alias{comet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize the EWAS results in a genomic region of interest
}
\description{
   The coMET is a R-based package to visualize the EWAS (epigenome-wide association scans) results in a genomic region of interest. The main feature of coMET is to plot the methylation correlation between CpG sites and demonstrate the significance level of EWAS results in the selected region. coMET package generates the plots of association, co-methylation patterns and a series of annotation tracks at genomic scale.
}
\usage{
comet(MYDATA.FILE = NULL, MYDATA.FORMAT = "SITE", MYDATA.LARGE.FILE = NULL, MYDATA.LARGE.FORMAT = "SITE", CORMATRIX.FILE = NULL, CORMATRIX.METHOD = "spearman", CORMATRIX.FORMAT = "RAW", CORMATRIX.COLOR.SCHEME = "bluewhitered", MYDATA.REF = NULL, START = NULL, END = NULL, ZOOM = FALSE, LAB.Y = "log", PVAL.THRESHOLD = 1e-05, DISP.PVAL.THRESHOLD = 1, DISP.ASSOCIATION = FALSE, DISP.ASSOCIATION.LARGE = FALSE, DISP.REGION = FALSE, DISP.REGION.LARGE = FALSE, SYMBOLS = "circle-fill", SYMBOLS.LARGE = NA, SAMPLE.LABELS = NULL, SAMPLE.LABELS.LARGE = NULL, USE.COLORS = TRUE , DISP.COLOR.REF = TRUE, COLOR.LIST = NULL, COLOR.LIST.LARGE = NULL, DISP.MYDATA = TRUE, BIOFEAT.USER.FILE = NULL, BIOFEAT.USER.TYPE = NULL, BIOFEAT.USER.TYPE.PLOT = NULL, GENOME = "hg19", DATASET.GENE = "hsapiens_gene_ensembl", DATASET.SNP = "hsapiens_snp", VERSION.DBSNP = "snp138", DATASET.SNP.STOMA = "hsapiens_snp_som", DATASET.REGULATION = "hsapiens_feature_set", DATASET.STRU = "hsapiens_structvar", DATASET.STRU.STOMA = "hsapiens_structvar_som", PATTERN.REGULATION = "GM12878", BROWSER.SESSION = "UCSC", TRACKS.GVIZ = NULL, TRACKS.GGBIO = NULL, TRACKS.TRACKVIEWER = NULL, DISP.MYDATA.NAMES = TRUE, DISP.COLOR.BAR = TRUE, DISP.PHYS.DIST = TRUE, DISP.LEGEND = TRUE, DISP.MARKER.LINES = TRUE, DISP.CORMATRIXMAP = TRUE, DISP.PVALUEPLOT =TRUE, DISP.TYPE = "symbol", DISP.MULT.LAB.X = FALSE, DISP.CONNECTING.LINES = TRUE, PALETTE.FILE = NULL, IMAGE.TITLE = NULL, IMAGE.NAME = "coMET", IMAGE.TYPE = NULL, IMAGE.SIZE = 3.5, FONT.FACTOR = NULL, SYMBOL.FACTOR = NULL, PRINT.IMAGE = TRUE, CONNECTING.LINES.FACTOR = 1.5, CONNECTING.LINES.ADJ = 0.01, CONNECTING.LINES.VERT.ADJ = -1, CONNECTING.LINES.FLEX = 0, config.file = NULL, VERBOSE = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MYDATA.FILE}{
 Name of a info file for coMET parameters
}
  \item{MYDATA.FORMAT}{
  Format of the input MYDATA.FILE. There are 6 different options: SITE, REGION, SITE_ASSOC, REGION_ASSOC
}
  \item{MYDATA.LARGE.FILE}{
  Name of additional info files for coMET parameters. Multiple files are acceptable, files should be separated by comma
}
  \item{MYDATA.LARGE.FORMAT}{
  Format of additional data visualised in p-value plot. Multiple files are acceptable, files should be separated by comma
}
  \item{CORMATRIX.FILE}{
 Name of raw data or pre-computed correlation matrix
}
  \item{CORMATRIX.METHOD}{
 Options for generating correlation matrix. There are three options: spearman, pearson and kendall
}
 \item{CORMATRIX.FORMAT}{
 Format of the input CORMATRIX.FILE. TThere are two options: raw files (RAW) and and pre-computed correlation matrix (CORMATRIX)
}
  \item{CORMATRIX.COLOR.SCHEME}{
  Color scheme options: heat, bluewhitered, cm, topo, gray, bluetored
}
  \item{MYDATA.REF}{
 The name of reference site listed in the MYDATA.FILE
}
  \item{START}{
  The first nucleotide position visualised
}
\item{END}{
   the last nucleotide position visualised
}
  \item{ZOOM}{
 Default=False
}
  \item{LAB.Y}{
 Scale of y-axis. Options: log or ln
}
 \item{PVAL.THRESHOLD}{
 Threshold of the significance. Displayed as a red dash line
}
  \item{DISP.PVAL.THRESHOLD}{
 Display only the findings pass PVAL.THRESHOLD
}
\item{DISP.ASSOCIATION}{
  Optional if MYDATA.FILE= SITE_ASSOC or REGION_ASSOC and contains the effect directions. A value which can be TRUE or False. If it is False (the default), the comethylation pattern is shown in p-value plot; if it is True, the effect directions is shown on the plot
}
  \item{DISP.ASSOCIATION.LARGE}{
 Optinal if MYDATA.LARGE.FILE contains the effect directions (SITE_ASSOC, REGION_ASSOC). A value which can be TRUE or False. If it is False (the default), the comethylation pattern is shown in p-value plot; if it is True, the effect directions is shown on the plot
}
  \item{DISP.REGION}{
  Optional if MYDATA.FILE= REGION or REGION_ASSOC. A value which can be TRUE or False (the default). If it is TRUE, the region of genomic element will be shown by a continuous line with the color of element, in addition to a symbole in center of region
}
  \item{DISP.REGION.LARGE}{
  Optional if MYDATA.LARGE.FILE= REGION or REGION_ASSOC. A value which can be TRUE or False (the default). If it is TRUE, the region of genomic element will be shown by a continuous line with the color of element, in addition to a symbole in center of region
}
  \item{SYMBOLS}{
  The symbol shown in the p-value plot. Options: circle, square, diamond, triangle. Symbols can be filled by appending -fill, e.g. square-fill. In this case, the SNP.FILE ASSOC column is read and an up-triangle and down-triangle are used to indicate positive and negative association (indicating susceptibility or protective alleles), respectively. Example: circle,NA,diamond-fill,triangle
}
  \item{SYMBOLS.LARGE}{
  The symbol to visualise the data defined in MYDATA.LARGE.FILE. Options: circle, square, diamond, triangle; Symbols can either be filled or not filled by appending -fill e.s., square-fill. NA may be specified. In this case, the SNP.FILE ASSOC column is read and an up-triangle and down-triangle are used to indicate positive and negative association (indicating susceptibility or protective alleles), respectively. Example: circle,NA,diamond-fill,triangle
}
  \item{SAMPLE.LABELS}{
  Labels for sample described in MYDATA.FILE to put in the legend
}
  \item{SAMPLE.LABELS.LARGE}{
  Labels for sample described in MYDATA.LARGE.FILE to put in the legend
}
  \item{USE.COLORS}{
  Use the colors defined or grey color
}
  \item{DISP.COLOR.REF}{
  True if you want to color connection line related to the reference probe in purple, FALSE if you do not want to color connection line related to the reference probe
}
  \item{COLOR.LIST}{
  List of colors for displaying p-value symbols related to data coming from MYDATA.FILE
}
  \item{COLOR.LIST.LARGE}{
  List of colors for displaying p-value symbols related to data coming from MYDATA.LARGE.FILE
}
  \item{DISP.MYDATA}{
  Logical True or False. The p-value plot is shown if it is TRUE; else the plot will be defined by GViz
}
 \item{BIOFEAT.USER.FILE}{
  Name of data file to visualise in tracks. They are separated by comma.
}
  \item{BIOFEAT.USER.TYPE}{
  Track type, multiple tracks can be shown (comma-separated): DataTrack, AnnotationTrack, GeneRegionTrack.
}
  \item{BIOFEAT.USER.TYPE.PLOT}{
  Format of plot if the data are DataTrack. They are separated by comma.
}
  \item{GENOME}{
  The human genome reference file. e.g. hg19 for Human genome 19 (NCBI 37)
}
  \item{DATASET.GENE}{
  The gene names from ENSEMBL. e.g. hsapiens_gene
}
  \item{DATASET.SNP}{
  Name of SNP database from ENSEMBL; Default: hsapiens_snp
}
  \item{VERSION.DBSNP}{
  Name of dbSNP used; Default : snp138 version from DBSNP
}
 \item{DATASET.SNP.STOMA}{
  Optional. Name of somatic SNP database from ENSEMBL. Default : hsapiens_snp_som
}
  \item{DATASET.REGULATION}{
  Optional. Name of regulation database from ENSEMBL. Default: hsapiens_feature_set
}
  \item{DATASET.STRU}{
  Optional. Name of structural variation database from ENSEMBL. Default: hsapiens_structvar
}
  \item{DATASET.STRU.STOMA}{
  Optional. Name of somatic structural variation database from ENSEMBL. Default: hsapiens_structvar_som
}
  \item{PATTERN.REGULATION}{
  Name of cell type of DATASET.REGULATION. Default: GM12878
}
  \item{BROWSER.SESSION}{
  Name of database for BioMART connection. Default: UCSC
}
 \item{TRACKS.GVIZ}{
  list of tracks created by Gviz.
}
  \item{TRACKS.GGBIO}{
  list of tracks created by ggbio.
}
  \item{TRACKS.TRACKVIEWER}{
  list of tracks created by track viewer.
}
 \item{DISP.MYDATA.NAMES}{
  Logical. If it is True (default), it displays the name of CpG sites.
}
  \item{DISP.COLOR.BAR}{
  color legend for the correlation matrix (range -1 to 1). Default: blue-white-red
}
  \item{DISP.PHYS.DIST}{
  display the distance of DNA sequence on plots instead of correlation matrix
}
\item{DISP.LEGEND}{
   Logical. Display the sample labels and corresponding symbols on the lower right side
}
  \item{DISP.MARKER.LINES}{
  Logical: True (default) or False. If it is False the red line for PVAL.THRESHOLD is not shown
}
  \item{DISP.CORMATRIXMAP}{
  Logical: True (default) or False. If it is False,  correlation matrix is not shown
}
  \item{DISP.PVALUEPLOT}{
  Logical: True (default) or False. If it is False,  pvalue plot is not shown
}
  \item{DISP.TYPE}{
  Default: symbol
}
  \item{DISP.MULT.LAB.X}{
  Logical. Display evenly spaced X-axis thick-labels; up to 5 labels are shown
}
  \item{DISP.CONNECTING.LINES}{
  Logical: True (default) or False. Display connecting lines between p-value plot and correlation matrix
}
  \item{PALETTE.FILE}{
  The path of file that contains color codes for heatmap. Colors are hexidecimal HTML color codes; one color per line; if do not want use this option, use the color defined by the option CORMATRIX.COLOR.SCHEME
}
  \item{IMAGE.TITLE}{
  Title of the plot
}
  \item{IMAGE.NAME}{
  Name of the plot file
}
  \item{IMAGE.TYPE}{
  Options: pdf or eps
}
  \item{IMAGE.SIZE}{
  Default:3.5 inches. Possible sizes : 3.5 or 7
}
  \item{FONT.FACTOR}{
  Font size of the sample labels. Range: 0-1
}
  \item{SYMBOL.FACTOR}{
  Size of the symbols. Range: 0-1
}
   \item{PRINT.IMAGE}{
  Print image in file or not.
}
  \item{CONNECTING.LINES.FACTOR}{
  Length of the connecting lines. Range: 0-2
}
   \item{CONNECTING.LINES.ADJ}{
  Position of the connecting lines horizontally. Negative values shift the connecting lines to the left and positive values shift the lines to the right. Range: (-1;1) option -1 means; no connecting lines adj
}
  \item{CONNECTING.LINES.VERT.ADJ}{
  Position of the connecting lines vertically. Can be used to vertically adjust the position of connecting lines in relation to cpg names. More negative value shift the connecting lines down. Range: -0.5 - 0, if -1,, use the default value related to the plot size (-0.5 for 3.5 plot size; -0.7 for 7.5 plot size)
}
  \item{CONNECTING.LINES.FLEX}{
   Adjusts the spread of the connecting lines. Range: 0-2
}
  \item{config.file}{
    Configuration file that contains the values of options instead of defining by command line 
}
   \item{VERBOSE}{
    DEFAULT=FALSE. If it is TRUE, it shows the comments.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Create a plot
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
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4comet.txt")

chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"

if(interactive()){
  cat("interactive")
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
  cat("Non interactive")
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
comet(config.file=configfile,TRACKS.GVIZ=listgviz, 
      VERBOSE=FALSE, PRINT.IMAGE=FALSE,DISP.PVALUEPLOT=FALSE)
}



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }% __ONLY ONE__ keyword per line
