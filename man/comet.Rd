\name{comet}
\alias{comet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize EWAS results in a genomic region of interest
}
\description{
   coMET is an R-based package to visualize EWAS (epigenome-wide association scans) results in a genomic region of interest. The main feature of coMET is to plot the the significance level of EWAS results in the selected region, along with correlation in DNA methylation values between CpG sites in the region. The coMET package generates plots of phenotype-association, co-methylation patterns, and a series of annotation tracks.
}
\usage{
comet(MYDATA.FILE = NULL, MYDATA.FORMAT = "SITE", MYDATA.LARGE.FILE = NULL, 
      MYDATA.LARGE.FORMAT = "SITE", CORMATRIX.FILE = NULL, CORMATRIX.METHOD = "spearman", 
      CORMATRIX.FORMAT = "RAW", CORMATRIX.COLOR.SCHEME = "bluewhitered", MYDATA.REF = NULL, 
      START = NULL, END = NULL, ZOOM = FALSE, LAB.Y = "log", PVAL.THRESHOLD = 1e-05, 
      DISP.PVAL.THRESHOLD = 1, DISP.ASSOCIATION = FALSE, DISP.ASSOCIATION.LARGE = FALSE, 
      DISP.REGION = FALSE, DISP.REGION.LARGE = FALSE, SYMBOLS = "circle-fill", 
      SYMBOLS.LARGE = NA, SAMPLE.LABELS = NULL, SAMPLE.LABELS.LARGE = NULL, 
      USE.COLORS = TRUE , DISP.COLOR.REF = TRUE, COLOR.LIST = NULL, COLOR.LIST.LARGE = NULL, 
      DISP.MYDATA = TRUE, BIOFEAT.USER.FILE = NULL, BIOFEAT.USER.TYPE = NULL, 
      BIOFEAT.USER.TYPE.PLOT = NULL, GENOME = "hg19", DATASET.GENE = "hsapiens_gene_ensembl", 
      TRACKS.GVIZ = NULL, TRACKS.GGBIO = NULL, TRACKS.TRACKVIEWER = NULL, 
      DISP.MYDATA.NAMES = TRUE, DISP.COLOR.BAR = TRUE, DISP.PHYS.DIST = TRUE, 
      DISP.LEGEND = TRUE, DISP.MARKER.LINES = TRUE, DISP.CORMATRIXMAP = TRUE, 
      DISP.PVALUEPLOT =TRUE, DISP.TYPE = "symbol", DISP.MULT.LAB.X = FALSE, 
      DISP.CONNECTING.LINES = TRUE, PALETTE.FILE = NULL, IMAGE.TITLE = NULL, 
      IMAGE.NAME = "coMET", IMAGE.TYPE = NULL, IMAGE.SIZE = 3.5, FONT.FACTOR = NULL, 
      SYMBOL.FACTOR = NULL, PRINT.IMAGE = TRUE, CONNECTING.LINES.FACTOR = 1.5, 
      CONNECTING.LINES.ADJ = 0.01, CONNECTING.LINES.VERT.ADJ = -1, 
      CONNECTING.LINES.FLEX = 0, config.file = NULL, VERBOSE = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MYDATA.FILE}{
 Name of the info file describing the coMET parameters
}
  \item{MYDATA.FORMAT}{
  Format of the input data in MYDATA.FILE. There are 4 different options: SITE, REGION, SITE_ASSOC, REGION_ASSOC
}
  \item{MYDATA.LARGE.FILE}{
  Name of additional info files describing the coMET parameters. File names should be comma-separated.
}
  \item{MYDATA.LARGE.FORMAT}{
  Format of additional data to be visualised in the p-value plot. File names should be comma-separated.
}
  \item{CORMATRIX.FILE}{
 Name of the raw data file or the pre-computed correlation matrix file.
}
  \item{CORMATRIX.METHOD}{
 Options for calculating the correlation matrix: spearman, pearson and kendall
}
 \item{CORMATRIX.FORMAT}{
 Format of the input CORMATRIX.FILE. TThere are two options: raw file (RAW if CpG sites are by column and samples by row or RAW_REV if CpG site are by row and samples by column) and pre-computed correlation matrix (CORMATRIX)
}
  \item{CORMATRIX.COLOR.SCHEME}{
  Color scheme options: heat, bluewhitered, cm, topo, gray, bluetored
}
  \item{MYDATA.REF}{
 The name of the reference CpG-site listed in MYDATA.FILE
}
  \item{START}{
  The first nucleotide position to be visualised
}
\item{END}{
   the last nucleotide position to be visualised
}
  \item{ZOOM}{
 Default=False
}
  \item{LAB.Y}{
 Scale of the y-axis. Options: log or ln
}
 \item{PVAL.THRESHOLD}{
 Significance threshold to be displayed as a red dashed line
}
  \item{DISP.PVAL.THRESHOLD}{
 Display only the findings that pass the value put in DISP.PVAL.THRESHOLD
}
\item{DISP.ASSOCIATION}{
  This logical option works only if MYDATA.FILE contains the effect direction (MYDATA.FORMAT=SITE_ASSOC or REGION_ASSOC). The value can be TRUE or FALSE: if FALSE (default), for each point of data in the p-value plot, the color of symbol is the color of co-methylation pattern between the point and the reference site; if TRUE, the effect direction is shown. If the association is positive, the color is the one defined with the option COLOR.LIST. On the other hand, if the association is negative, the color is the opposed color. 
}
  \item{DISP.ASSOCIATION.LARGE}{
  This logical option works only if MYDATA.LARGE.FILE contains the effect direction (MYDATA.LARGE.FORMAT=SITE_ASSOC or REGION_ASSOC). The value can be TRUE or FALSE: if FALSE (default), for each point of data in the p-value plot, the color of symbol is the color of co-methylation pattern between the point and the reference site; if TRUE, the effect direction is shown. If the association is positive, the color is the one defined with the option COLOR.LIST.LARGE. On the other hand, if the association is negative, the color is the opposed color.
}
  \item{DISP.REGION}{
  This logical option works only if MYDATA.FILE contains regions (MYDATA.FORMAT=REGION or REGION_ASSOC). The value can be TRUE or FALSE (default). If TRUE, the genomic element will be shown by a continuous line with the color of the element, in addition to the symbol at the center of the region. If FALSE, only the symbol is shown.
}
  \item{DISP.REGION.LARGE}{
  This logical option works only if MYDATA.LARGE.FILE contains regions (MYDATA.LARGE.FORMAT=REGION or REGION_ASSOC). The value can be TRUE or FALSE (default). If TRUE, the genomic element will be shown by a continuous line with the color of the element, in addition to the symbol at the center of the region. If FALSE, only the symbol is shown.
}
  \item{SYMBOLS}{
  The symbol shown in the p-value plot. Options: circle, square, diamond, triangle. Symbols can be filled by appending -fill, e.g. square-fill. Example: circle,diamond-fill,triangle
}
  \item{SYMBOLS.LARGE}{
  The symbol to visualise the data defined in MYDATA.LARGE.FILE. Options: circle, square, diamond, triangle; Symbols can either be filled or not filled by appending -fill e.s., square-fill. Example: circle,diamond-fill,triangle
}
  \item{SAMPLE.LABELS}{
  Labels for the sample described in MYDATA.FILE to include in the legend
}
  \item{SAMPLE.LABELS.LARGE}{
  Labels for the sample described in MYDATA.LARGE.FILE to include in the legend
}
  \item{USE.COLORS}{
  Use the colors defined or use the grey color scheme
}
  \item{DISP.COLOR.REF}{
  Logical option TRUE or FALSE (TRUE default). if TRUE, the connection line related to the reference probe is in purple, if FALSE if the connection line related to the reference probe stay black.
}
  \item{COLOR.LIST}{
  List of colors for displaying the P-value symbols related to the data in MYDATA.FILE
}
  \item{COLOR.LIST.LARGE}{
  List of colors for displaying the P-value symbols related to the data in MYDATA.LARGE.FILE
}
  \item{DISP.MYDATA}{
  logical option TRUE or FALSE. TRUE (default). If TRUE, the P-value plot is shown; if FALSE the plot will be defined by GViz
}
 \item{BIOFEAT.USER.FILE}{
  Name of data file to visualise in the tracks. File names should be comma-separated.
}
  \item{BIOFEAT.USER.TYPE}{
  Track type, where multiple tracks can be shown (comma-separated): DataTrack, AnnotationTrack, GeneRegionTrack.
}
  \item{BIOFEAT.USER.TYPE.PLOT}{
  Format of the plot if the data are shown with the Gviz's function called DataTrack (comma-separated)
}
  \item{GENOME}{
  The human genome reference file. e.g. "hg19" for Human genome 19 (NCBI 37), "grch37" (GRCh37),"grch38" (GRCh38)
}
  \item{DATASET.GENE}{
  The gene names from ENSEMBL. e.g. hsapiens_gene
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
  logical option TRUE or FALSE. If True (default), the names of the CpG sites are displayed.
}
  \item{DISP.COLOR.BAR}{
  Color legend for the correlation matrix (range -1 to 1). Default: blue-white-red
}
  \item{DISP.PHYS.DIST}{
  logical option (TRUE or FALSE). TRUE (default).Display the bp distance on the plots
}
\item{DISP.LEGEND}{
  logical option TRUE or FALSE. TRUE (default) Display the sample labels and corresponding symbols on the lower right side
}
  \item{DISP.MARKER.LINES}{
  logical option TRUE or FALSE. TRUE (default), if FALSE the red line for PVAL.THRESHOLD is not shown
}
  \item{DISP.CORMATRIXMAP}{
  logical option TRUE or FALSE. TRUE (default), if FALSE correlation matrix is not shown
}
  \item{DISP.PVALUEPLOT}{
  logical option (TRUE or FALSE). TRUE (default), if FALSE the pvalue plot is not shown
}
  \item{DISP.TYPE}{
  Default: symbol
}
  \item{DISP.MULT.LAB.X}{
  logical option TRUE or FALSE. FALSE (default).Display evenly spaced X-axis labels; up to 5 labels are shown.
}
  \item{DISP.CONNECTING.LINES}{
  logical option TRUE or FALSE. TRUE (default) displays connecting lines between p-value plot and correlation matrix
}
  \item{PALETTE.FILE}{
  File that contains color scheme for the heatmap. Colors are hexidecimal HTML color codes; one color per line; if you do not want to use this option, use the color defined by the option CORMATRIX.COLOR.SCHEME
}
  \item{IMAGE.TITLE}{
  Title of the plot
}
  \item{IMAGE.NAME}{
  Path and Name of the plot file without extension
}
  \item{IMAGE.TYPE}{
  Options: pdf or eps
}
  \item{IMAGE.SIZE}{
  Default: 3.5 inches. Possible sizes : 3.5 or 7
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
  Position of the connecting lines horizontally. Negative values shift the connecting lines to the left and positive values shift the lines to the right. Range: (-1;1) option -1 means no connecting lines.
}
  \item{CONNECTING.LINES.VERT.ADJ}{
  Position of the connecting lines vertically. Can be used to vertically adjust the position of the connecting lines in relation to the CpG-site names. Negative value shift the connecting lines down. Range: (-0.5 - 0), option -1 mean the default value related to the plot size (-0.5 for 3.5 plot size; -0.7 for 7.5 plot size)
}
  \item{CONNECTING.LINES.FLEX}{
   Adjusts the spread of the connecting lines. Range: 0-2
}
  \item{config.file}{
    Configuration file that contains the values of these options instead of defining these by command line 
}
   \item{VERBOSE}{
    logical option TRUE or FALSE. TRUE (default). If TRUE, shows comments.
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
myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")

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
comet(config.file=configfile, MYDATA.FILE=myinfofile, CORMATRIX.FILE=mycorrelation, 
      MYDATA.LARGE.FILE=myexpressfile, TRACKS.GVIZ=listgviz, 
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
comet(config.file=configfile, MYDATA.FILE=myinfofile, CORMATRIX.FILE=mycorrelation, 
      MYDATA.LARGE.FILE=myexpressfile, TRACKS.GVIZ=listgviz, 
      VERBOSE=FALSE, PRINT.IMAGE=FALSE,DISP.PVALUEPLOT=FALSE)
}



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }% __ONLY ONE__ keyword per line
