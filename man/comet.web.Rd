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
comet.web(mydata.file = NULL, mydata.format = c("site", "region", "site_asso", "region_asso"), 
        mydata.large.file = NULL, 
        mydata.large.format = c("site", "region", "site_asso", "region_asso"), 
        cormatrix.file = NULL, cormatrix.method = c("spearman", "pearson", "kendall"), 
        cormatrix.format = c("cormatrix", "raw","raw_rev"), 
        cormatrix.color.scheme = "heat", cormatrix.conf.level=0.05, 
        cormatrix.sig.level= 1, cormatrix.adjust="none",mydata.ref = NULL,
        genome="hg19", start = NULL, end = NULL, zoom = FALSE, lab.Y = "log", 
        pval.threshold = 1e-07, disp.pval.threshold = 1, 
        disp.association= FALSE, disp.association.large = FALSE, 
        disp.region = FALSE, disp.region.large = FALSE, symbols = "circle-fill", 
        symbols.large = NA, sample.labels = NULL, sample.labels.large = NULL, 
        use.colors = TRUE, disp.color.ref = TRUE, color.list = NULL, 
        color.list.large = NULL, biofeat.user.file = NULL,
        biofeat.user.type = c("GeneRegion", "Annotation", "Data"), 
        biofeat.user.type.plot = NULL, 
        list.tracks = "geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP", 
        pattern.regulation = "GM12878",
        image.title = NULL, image.name = "coMET", image.type = c("pdf", "eps"), 
        image.size = 3.5, print.image = FALSE, config.file = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{.
    \item{mydata.file}{
 Name of the info file describing the coMET parameters. It is mandatory and has to be a file in tabular format with a header. Info file can be a list of CpG sites with/without Beta value (DNA methylation level) or direction sign. If it is a site file then it is mandatory to have the 4 columns as shown below with headers in the same order. Beta can be the 5th column(optional) and it can be either a numeric value (positive or negative values) or only direction sign ("+", "-"). The number of columns and their types are defined but the option mydata.format.
}
    \item{mydata.format}{
    Format of the input data in mydata.file. There are 4 different options: site, region, site_asso, region_asso.
}
    \item{mydata.large.file}{
    Name of additional info files describing the coMET parameters. File names should be comma-separated. It is optional, but if you add some, they need to be file(s) in tabular format with a header. Additional info file can be a list of CpG sites with/without Beta value (DNA methylation level) or direction sign. If it is a site file then it is mandatory to have the 4 columns as shown below with headers in the same order. Beta can be the 5th column(optional) and it can be either a numeric value (positive or negative values) or only direction sign ("+", "-"). The number of columns and their types are defined but the option mydata.large.format.
}
    \item{mydata.large.format}{
    Format of additional data to be visualised in the p-value plot. Format should be comma-separated.There are 4 different options for each file: site, region, site_asso, region_asso.
}
    \item{cormatrix.file}{
 Name of the raw data file or the pre-computed correlation matrix file. It is mandatory and has to be a file in tabular format with an header.
}
    \item{cormatrix.method}{
    A character string indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman", can be abbreviated.
}
    \item{cormatrix.format}{
    A character string indicating which format of the input cormatrix.file is to be used. There are three options: raw file (raw if CpG sites are by column and samples by row or row_rev if CpG site are by row and samples by column) and pre-computed correlation matrix (cormatrix)
}
    \item{cormatrix.color.scheme}{
    A character string indicating which Color scheme options is to be used: heat, bluewhitered, cm, topo, gray, bluetored
}
    \item{cormatrix.conf.level}{
    Confidence level for the returned confidence interval. Currently only used for the Pearson product moment correlation coefficient if there are at least 4 complete pairs of observations. Default value= 0.05
}
    \item{cormatrix.sig.level}{
    Significant level to visualise the correlation. If the correlation has a pvalue under the significant level, the correlation will be colored in "goshwhite", else the color is related to the correlation level and the color scheme choosen.Default value =1.
}
    \item{cormatrix.adjust}{
    indicates which adjustment for multiple tests should be used. "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".Default value="none"
}
 
    \item{mydata.ref}{
    The name of the reference omic feature (e.g. CpG-site) listed in mydata.file
}
    \item{genome}{
    The human genome reference file. e.g. "hg19" for Human genome 19 (NCBI 37), "grch37" (GRCh37),"grch38" (GRCh38)
}
    \item{start}{
    The first nucleotide position to be visualised. It could be bigger or smaller than the first position of our list of omic features.
}
    \item{end}{
    the last nucleotide position to be visualised. It has to be bigger than the value in the option start, but it could be smaller or bigger than the last position of our list of omic features.
}
    \item{zoom}{
    logical option TRUE or FALSE. FALSE (default)
}
    \item{lab.Y}{
    Scale of the y-axis. Options: log or ln
}
    \item{pval.threshold}{
    Significance threshold to be displayed as a red dashed line
}
    \item{disp.pval.threshold}{
    Display only the findings that pass the value put in disp.pval.threshold
}
    \item{disp.association}{
    This logical option works only if mydata.file contains the effect direction (mydata.format=site_asso or region_asso). The value can be TRUE or FALSE: if FALSE (default), for each point of data in the p-value plot, the color of symbol is the color of co-methylation pattern between the point and the reference site; if TRUE, the effect direction is shown. If the association is positive, the color is the one defined with the option color.list. On the other hand, if the association is negative, the color is the opposed color. 
}
    \item{disp.association.large}{
    This logical option works only if mydata.large.file contains the effect direction (MYDATA.large.FORMA=site_asso or region_asso). The value can be TRUE or FALSE: if FALSE (default), for each point of data in the p-value plot, the color of symbol is the color of co-methylation pattern between the point and the reference site; if TRUE, the effect direction is shown. If the association is positive, the color is the one defined with the option color.list.large. On the other hand, if the association is negative, the color is the opposed color.
}
    \item{disp.region}{
    This logical option works only if mydata.file contains regions (mydata.format=region or region_asso). The value can be TRUE or FALSE (default). If TRUE, the genomic element will be shown by a continuous line with the color of the element, in addition to the symbol at the center of the region. If FALSE, only the symbol is shown.
}
    \item{disp.region.large}{
    This logical option works only if mydata.large.file contains regions (mydata.large.format=region or region_asso). The value can be TRUE or FALSE (default). If TRUE, the genomic element will be shown by a continuous line with the color of the element, in addition to the symbol at the center of the region. If FALSE, only the symbol is shown.
}
    \item{symbols}{
    The symbol shown in the p-value plot. Options: circle, square, diamond, triangle. symbols can be filled by appending -fill, e.g. square-fill. Example: circle,diamond-fill,triangle
}
    \item{symbols.large}{
    The symbol to visualise the data defined in mydata.large.file. Options: circle, square, diamond, triangle; symbols can either be filled or not filled by appending -fill e.s., square-fill. Example: circle,diamond-fill,triangle
}
    \item{sample.labels}{
    Labels for the sample described in mydata.file to include in the legend
}
    \item{sample.labels.large}{
    Labels for the sample described in mydata.large.file to include in the legend
}
    \item{use.colors}{
    Use the colors defined or use the grey color scheme
}
    \item{disp.color.ref}{
    Logical option TRUE or FALSE (TRUE default). if TRUE, the connection line related to the reference probe is in purple, if FALSE if the connection line related to the reference probe stay black.
}
    \item{color.list}{
    List of colors for displaying the P-value symbols related to the data in mydata.file
}
    \item{color.list.large}{
    List of colors for displaying the P-value symbols related to the data in mydata.large.file
}
    \item{biofeat.user.file}{
    Name of data file to visualise in the tracks. File names should be comma-separated.
}
    \item{biofeat.user.type}{
    Track type, where multiple tracks can be shown (comma-separated): DataTrack, AnnotationTrack, GeneRegionTrack.
}
    \item{biofeat.user.type.plot}{
    Format of the plot if the data are shown with the Gviz's function called DataTrack (comma-separated)
}
    \item{list.tracks}{
    List of annotation tracks to visualise. Options include geneENSEMBL, CGI, ChromHMM, DNAse, RegENSEMBL, SNP, transcriptENSEMBL, SNPstoma, SNPstru, SNPstrustoma, ISCA, COSMIC, GAD, ClinVar, GeneReviews, GWAS, ClinVarCNV, GCcontent, genesUCSC, xenogenesUCSC. 
}
    \item{pattern.regulation}{
    The cell/tissue or the list of cells/tissues to visualise in the regulation region defined by Broad ChromHMM
}
    \item{image.title}{
    Title of the plot
}
    \item{image.name}{
    The path and the name of the plot file without extension. The extension will be added by coMET depending on the option image.type.
}
    \item{image.type}{
    Options: pdf or eps
}
    \item{image.size}{
    Default: 3.5 inches. Possible sizes : 3.5 or 7
}
    \item{print.image}{
    Print image in file or not.
}
    \item{config.file}{
    Configuration file contains the values of these options instead of defining these by command line. It is a file where each line is one option. The name of option is in capital and is separated to its value by "=". If there are multiple values such as for the option list.tracks or the options for additional data, you need to separated them by a "comma".
}
    \item{verbose}{
    logical option TRUE or FALSE. TRUE (default). If TRUE, shows comments.
}
}
\details{
The function is limited to visualize 120 omic features.
}
\value{
Create a plot in pdf or eps format depending to some options
%%    ~Describe the value returned
%%    If it is a LIST, use
%%    \item{comp1 }{Description of 'comp1'}
%%    \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
http://epigen.kcl.ac.uk/comet/
}
\author{
Tiphaine Martin
}
%%\note{
%%    ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{comet}},\code{\link{comet.list}}
}
\examples{
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4webserver.txt")
myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")

comet.web(config.file=configfile, mydata.file=myinfofile, cormatrix.file=mycorrelation, 
    mydata.large.file=myexpressfile, print.image=FALSE,verbose=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
