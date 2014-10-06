\name{comet.web}
\alias{comet.web}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize the EWAS results in a genomic region of interest with a
predefined annotation tracks
}
\description{
The coMET is a R-based package to visualize the EWAS (epigenome-wide association scans) results in a genomic region of interest. The main feature of coMET is to plot the methylation correlation between CpG sites and demonstrate the significance level of EWAS results in the selected region. coMET package generates the plots of association, co-methylation patterns and a series of annotation tracks at genomic scale.
}
\usage{
comet.web(MYDATA.FILE = NULL, MYDATA.FORMAT = c("SITE", "REGION", "SITE_ASSO", "REGION_ASSO"), MYDATA.LARGE.FILE = NULL, MYDATA.LARGE.FORMAT = c("SITE", "REGION", "SITE_ASSO", "REGION_ASSO"), CORMATRIX.FILE = NULL, CORMATRIX.METHOD = c("spearman", "pearson", "kendall"), CORMATRIX.FORMAT = c("CORMATRIX", "RAW"), CORMATRIX.COLOR.SCHEME = "heat", MYDATA.REF = NULL, START = NULL, END = NULL, ZOOM = FALSE, LAB.Y = "log", PVAL.THRESHOLD = 1e-07, DISP.PVAL.THRESHOLD = 1, DISP.ASSOCIATION = FALSE, DISP.ASSOCIATION.LARGE = FALSE, DISP.REGION = FALSE, DISP.REGION.LARGE = FALSE, SYMBOLS = "circle-fill", SYMBOLS.LARGE = NA, SAMPLE.LABELS = NULL, SAMPLE.LABELS.LARGE = NULL, USE.COLORS = TRUE, DISP.COLOR.REF = TRUE, COLOR.LIST = NULL, COLOR.LIST.LARGE = NULL, BIOFEAT.USER.FILE = NULL, BIOFEAT.USER.TYPE = c("GeneRegion", "Annotation", "Data"), BIOFEAT.USER.TYPE.PLOT = NULL, LIST.TRACKS = "geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP", PATTERN.REGULATION = "GM12878", IMAGE.TITLE = NULL, IMAGE.NAME = "coMET", IMAGE.TYPE = c("pdf", "eps"), IMAGE.SIZE = 3.5, PRINT.IMAGE = FALSE, config.file = NULL, VERBOSE = FALSE)
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
 File with correlation matrix or raw data to create the correlation
}
  \item{CORMATRIX.METHOD}{
 Options for generating correlation matrix. There are three options: spearman, pearson and kendall
}
  \item{CORMATRIX.FORMAT}{
 Format of the input CORMATRIX.FILE. Two options for raw files (RAW or DTR_RAW) and and pre-computed correlation matrix (CORMATRIX)
}
  \item{CORMATRIX.COLOR.SCHEME}{
  Color scheme options: heat, ETC
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
  Optional if MYDATA.FILE= REGION or REGION_ASSOC. A value which can be TRUE or False (the default). If it is TRUE, the plot will shown in the newly defined x-scale limits
}
  \item{DISP.REGION.LARGE}{
  Optional if MYDATA.LARGE.FILE= REGION or REGION_ASSOC. A value which
  can be TRUE or False (the default). If it is TRUE, the plot will shown
  in the newly defined x-scale limits
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
  DEFAULT= TRUE; it is FALSE, no color
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
   \item{BIOFEAT.USER.FILE}{
  Name of data file to visualise in tracks. They are separated by comma.
}
  \item{BIOFEAT.USER.TYPE}{
  Track type, multiple tracks can be shown (comma-separated): DataTrack, AnnotationTrack, GeneRegionTrack.
}
  \item{BIOFEAT.USER.TYPE.PLOT}{
  Format of plot if the data are DataTrack. They are separated by comma.
}
  \item{LIST.TRACKS}{
  Options from geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP
}
  \item{PATTERN.REGULATION}{
  The tissue or the list of tissues to visualise the regulation region defined by Broad ChromHMM
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
  \item{PRINT.IMAGE}{
DEFAULT=FALSE. if the value is false, it does not produce the plot in
  a file. If the value is true, it print the plot in a file
}
  \item{config.file}{
  Configuration file that contains the values of options instead of defining by command line 
}
  \item{VERBOSE}{
  Visualisation of comment to understand what happens
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
configfile <- file.path(extdata, "config_cyp1b1_zoom.txt")
comet.web(config.file=configfile,PRINT.IMAGE=FALSE,VERBOSE=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
