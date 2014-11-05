\name{comet.web}
\alias{comet.web}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize EWAS results in a genomic region of interest with predefined annotation tracks
}
\description{
   coMET is an R-based package to visualize EWAS (epigenome-wide association scans) results in a genomic region of interest. The main feature of coMET is to plot the the significance level of EWAS results in the selected region, along with correlation in DNA methylation values between CpG sites in the region. The coMET package generates plots of phenotype-association, co-methylation patterns, and a series of annotation tracks.
}
\usage{
comet.web(MYDATA.FILE = NULL, MYDATA.FORMAT = c("SITE", "REGION", "SITE_ASSO", "REGION_ASSO"), 
          MYDATA.LARGE.FILE = NULL, 
          MYDATA.LARGE.FORMAT = c("SITE", "REGION", "SITE_ASSO", "REGION_ASSO"), 
          CORMATRIX.FILE = NULL, CORMATRIX.METHOD = c("spearman", "pearson", "kendall"), 
          CORMATRIX.FORMAT = c("CORMATRIX", "RAW","RAW_REV"), 
          CORMATRIX.COLOR.SCHEME = "heat", MYDATA.REF = NULL,
          GENOME="hg19", START = NULL, END = NULL, ZOOM = FALSE, LAB.Y = "log", 
          PVAL.THRESHOLD = 1e-07, DISP.PVAL.THRESHOLD = 1, 
          DISP.ASSOCIATION = FALSE, DISP.ASSOCIATION.LARGE = FALSE, 
          DISP.REGION = FALSE, DISP.REGION.LARGE = FALSE, SYMBOLS = "circle-fill", 
          SYMBOLS.LARGE = NA, SAMPLE.LABELS = NULL, SAMPLE.LABELS.LARGE = NULL, 
          USE.COLORS = TRUE, DISP.COLOR.REF = TRUE, COLOR.LIST = NULL, 
          COLOR.LIST.LARGE = NULL, BIOFEAT.USER.FILE = NULL,
          BIOFEAT.USER.TYPE = c("GeneRegion", "Annotation", "Data"), 
          BIOFEAT.USER.TYPE.PLOT = NULL, 
          LIST.TRACKS = "geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP", 
          PATTERN.REGULATION = "GM12878",
          IMAGE.TITLE = NULL, IMAGE.NAME = "coMET", IMAGE.TYPE = c("pdf", "eps"), 
          IMAGE.SIZE = 3.5, PRINT.IMAGE = FALSE, config.file = NULL, VERBOSE = FALSE)
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
 Format of the input CORMATRIX.FILE. There are three options: raw file (RAW if CpG sites are by column and samples by row or RAW_REV if CpG site are by row and samples by column) and pre-computed correlation matrix (CORMATRIX)
}
  \item{CORMATRIX.COLOR.SCHEME}{
  Color scheme options: heat, bluewhitered, cm, topo, gray, bluetored
}
  \item{MYDATA.REF}{
 The name of the reference CpG-site listed in MYDATA.FILE
}
  \item{GENOME}{
  The human genome reference file. e.g. "hg19" for Human genome 19 (NCBI 37), "grch37" (GRCh37),"grch38" (GRCh38)
}
  \item{START}{
  The first nucleotide position to be visualised
}
\item{END}{
   the last nucleotide position to be visualised
}
  \item{ZOOM}{
 logical option TRUE or FALSE. FALSE (default)
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
   This logical option works only if MYDATA.LARGE.FILE contains the effect direction (MYDATA.LARGE.FORMA=SITE_ASSOC or REGION_ASSOC). The value can be TRUE or FALSE: if FALSE (default), for each point of data in the p-value plot, the color of symbol is the color of co-methylation pattern between the point and the reference site; if TRUE, the effect direction is shown. If the association is positive, the color is the one defined with the option COLOR.LIST.LARGE. On the other hand, if the association is negative, the color is the opposed color.
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
 \item{BIOFEAT.USER.FILE}{
  Name of data file to visualise in the tracks. File names should be comma-separated.
}
  \item{BIOFEAT.USER.TYPE}{
  Track type, where multiple tracks can be shown (comma-separated): DataTrack, AnnotationTrack, GeneRegionTrack.
}
  \item{BIOFEAT.USER.TYPE.PLOT}{
  Format of the plot if the data are shown with the Gviz's function called DataTrack (comma-separated)
}
  \item{LIST.TRACKS}{
   List of annotation tracks to visualise. Options include geneENSEMBL, CGI, ChromHMM, DNAse, RegENSEMBL, SNP, transcriptENSEMBL, SNPstoma, SNPstru, SNPstrustoma, ISCA, COSMIC, GAD, ClinVar, GeneReviews, GWAS, ClinVarCNV, GCcontent, genesUCSC, xenogenesUCSC. 
}
  \item{PATTERN.REGULATION}{
  The cell/tissue or the list of cells/tissues to visualise in the regulation region defined by Broad ChromHMM
}
  \item{IMAGE.TITLE}{
  Title of the plot
}
  \item{IMAGE.NAME}{
  path and Name of the plot file without extension
}
  \item{IMAGE.TYPE}{
  Options: pdf or eps
}
  \item{IMAGE.SIZE}{
  Default: 3.5 inches. Possible sizes : 3.5 or 7
}
   \item{PRINT.IMAGE}{
  Print image in file or not.
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
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4webserver.txt")
myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")

comet.web(config.file=configfile,  MYDATA.FILE=myinfofile, CORMATRIX.FILE=mycorrelation, 
      MYDATA.LARGE.FILE=myexpressfile, PRINT.IMAGE=FALSE,VERBOSE=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
